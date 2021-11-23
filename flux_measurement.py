import numpy as np
import sys
from astropy.io import fits
from scipy.optimize import curve_fit
from photutils.detection import find_peaks
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, SourceCatalog, make_source_mask
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.background import Background2D
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
import pandas as pd


def gaussian(xy, peak, Fwhm, offset, x0, y0):
    x, y = xy
    sigma = Fwhm / (2 * np.sqrt(2 * np.log(2)))
    return peak * np.exp(-((y - y0) ** 2 + (x - x0) ** 2) / (2 * sigma ** 2)) + offset


def distance2(loc1, loc2):
    x0, y0, _ = loc1
    x1, y1, _ = loc2
    return (x1 - x0) ** 2 + (y1 - y0) ** 2


# INITIAL SET UP #
if len(sys.argv) == 4:
    _, img_file, output_file, searchradius = sys.argv
    searchradius = float(searchradius)
    try:
        assert int(searchradius) == searchradius
    except AssertionError:
        raise ValueError('Search radius needs to be an integer')
    searchradius = int(searchradius)
elif len(sys.argv) == 3:
    _, img_file, output_file = sys.argv
    searchradius = 10
else:
    raise ValueError('Incorrect number of arguments specified. First argument should be the filename for the image '
                     'and the second argument should be the output file.')

with fits.open(img_file) as f:
    data = f[0].data  # 1024 x 1024
data[np.where(np.isnan(data))] = 0
mask = make_source_mask(data, nsigma=15, npixels=5, dilate_size=11)
bckg = Background2D(data=data, box_size=4, coverage_mask=mask, fill_value=0.0).background
data = data - bckg
threshold = detect_threshold(data=data, nsigma=15, mask_value=0.0)
# END OF SET UP #

# GOAL OF THIS SECTION IS TO MEASURE FWHM OF THE IMAGE #
locs_table = find_peaks(data=data, threshold=threshold, npeaks=100)  # brightest 100 spots in image

assert locs_table is not None

locs = [(xloc, yloc, peak) for xloc, yloc, peak in
        zip(locs_table['x_peak'], locs_table['y_peak'], locs_table['peak_value'])]
estimated_saturated_point = np.percentile(locs_table['peak_value'], 95)

fwhm_measurements = list()

for _, loc in enumerate(locs):  # might use an index later on to modify loc before next round if doing across images
    x0, y0, estimated_peak = loc
    if np.min([distance2(loc1, loc) for loc1 in locs if loc1 != loc]) < 400:  # within 20 pixel radius
        continue
    if estimated_peak > estimated_saturated_point:
        continue
    searchbox = data[y0 - searchradius: y0 + searchradius + 1, x0 - searchradius: x0 + searchradius + 1]

    # building out arrays so that we get a full coordinate system when they are zipped together
    uppery, upperx = searchbox.shape
    yvals, xvals = np.arange(0, uppery, 1.0) - searchradius, np.arange(0, upperx, 1.0) - searchradius
    numyvals, numxvals = len(yvals), len(xvals)
    yvals = np.array(list(yvals) * numxvals)
    xvals = np.array([[xval] * numyvals for xval in xvals]).flatten()

    vals = [searchbox[int(y + searchradius), int(x + searchradius)] for y, x in zip(yvals, xvals)]
    vals_w_distances = pd.DataFrame({'vals': vals, 'distance squared': [y ** 2 + x ** 2 for y, x in zip(yvals, xvals)],
                                     'y': yvals, 'x': xvals})

    # estimating FWHM
    avg_vals_at_distance = {d2: np.mean(vals_w_distances[vals_w_distances['distance squared'] == d2]['vals'])
                            for d2 in np.unique(vals_w_distances['distance squared'])}
    below_half = [d2 for d2 in avg_vals_at_distance.keys() if avg_vals_at_distance[d2] <= 0.5 * estimated_peak]
    above_half = [d2 for d2 in avg_vals_at_distance.keys() if avg_vals_at_distance[d2] >= 0.5 * estimated_peak]
    if len(set(below_half).intersection(above_half)) > 0:
        fd2 = np.mean(list(set(below_half).intersection(above_half)))
    else:
        fd2 = np.mean([np.min(below_half), np.max(above_half)])
    estimated_fwhm = 2 * np.sqrt(fd2)

    # background has been subtracted so offset should be close to zero
    estimated_offset = 0

    # getting new estimate on FWHM
    finalsearchradius = int(np.ceil(estimated_fwhm))
    finalsearchbox = data[y0 - finalsearchradius: y0 + finalsearchradius + 1,
                          x0 - finalsearchradius: x0 + finalsearchradius + 1]
    uppery, upperx = finalsearchbox.shape
    yvals, xvals = np.arange(0, uppery, 1.0) - finalsearchradius, np.arange(0, upperx, 1.0) - finalsearchradius
    numyvals, numxvals = len(yvals), len(xvals)
    yvals = np.array(list(yvals) * numxvals)
    xvals = np.array([[xval] * numyvals for xval in xvals]).flatten()
    try:
        finalvals = [searchbox[int(y + finalsearchradius), int(x + finalsearchradius)] for y, x in zip(yvals, xvals)]
    except IndexError:
        continue

    guesses = [estimated_peak * 0.9, estimated_fwhm, estimated_offset, 0, 0]
    bounds = ((0, 0, 0, -1 * estimated_fwhm / 2, -1 * estimated_fwhm / 2),
              (np.max(finalsearchbox) * 1.1, np.inf, np.max(finalsearchbox), estimated_fwhm / 2, estimated_fwhm / 2))
    try:
        optimalparams, covariance_matrix = curve_fit(f=gaussian, xdata=(yvals, xvals), ydata=finalvals,
                                                     p0=guesses, bounds=bounds)
    except RuntimeError:  # if SciPy can't find params and gives up
        continue
    optimalpeak, optimalfwhm, optimaloffset, optimalx0, optimaly0 = optimalparams
    fwhm_measurements.append(optimalfwhm)

# average of measurements excluding outliers
FWHM = np.mean([fm for fm in fwhm_measurements
                if np.percentile(fwhm_measurements, 3) < fm < np.percentile(fwhm_measurements, 97)])
# END OF FINDING IMAGE FWHM #

# GOAL OF THIS SECTION IS TO ACTUALLY MEASURE STAR FLUXES #
# get positions (we aren't using the ones from before since that was a more limited subset of stars)
sigma = FWHM * gaussian_fwhm_to_sigma
kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
kernel.normalize()
srcs = detect_sources(data, threshold, npixels=5, kernel=kernel)
dsrcs = deblend_sources(data, srcs, npixels=5, kernel=kernel)
cat = SourceCatalog(data, dsrcs)
cat_table = cat.to_table()
positions = [(x, y) for x, y in zip(cat_table['xcentroid'], cat_table['ycentroid'])]

# measure fluxes (normalized by area)
aperature = CircularAperture(positions, r=3*FWHM)
flux_table = aperture_photometry(data, aperature)
xcenters, ycenters, fluxes = flux_table['xcenter'], flux_table['ycenter'], flux_table['aperture_sum']
fluxes = np.array(fluxes) / (np.pi * ((3 * FWHM) ** 2))  # normalizing by area
dataframe_data = {'X': xcenters, 'Y': ycenters, 'Flux': fluxes}

df = pd.DataFrame(data=dataframe_data)
df.to_csv(output_file, index=False)
