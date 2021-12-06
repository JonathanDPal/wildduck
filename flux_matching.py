import numpy as np
import pandas as pd
import astroalign as aa
from astropy.io import fits
import sys

# USER INPUTS #
_, direc, fwhm = sys.argv
# Original images that flux measurements were made on
r_imagefile = f'{direc}/{direc}-R-stack.fit'
b_imagefile = f'{direc}/{direc}-B-stack.fit'
v_imagefile = f'{direc}/{direc}-V-stack.fit'
# CSV files that were output of 'flux_measurement.py'
r_fluxtable = f'{direc}/{direc}-Rfluxes.csv'
b_fluxtable = f'{direc}/{direc}-Bfluxes.csv'
v_fluxtable = f'{direc}/{direc}-Vfluxes.csv'
FWHM = float(fwhm)
num_sat_stars = 5  # CANNOT BE ZERO
output_file = f'{direc}/{direc}-fluxes.csv'
# END OF USER INPUTS #


def distance2(loc1, loc2):
    """
    Distance squared
    """
    x0, y0 = loc1
    x1, y1 = loc2
    return (x1 - x0) ** 2 + (y1 - y0) ** 2


v_fluxes = pd.read_csv(v_fluxtable).sort_values('Flux', ascending=True)[:-1 * num_sat_stars]
r_fluxes = pd.read_csv(r_fluxtable).sort_values('Flux', ascending=True)[:-1 * num_sat_stars]
b_fluxes = pd.read_csv(b_fluxtable).sort_values('Flux', ascending=True)[:-1 * num_sat_stars]

with fits.open(r_imagefile) as f:
    r_img = f[0].data
with fits.open(b_imagefile) as f:
    b_img = f[0].data
with fits.open(v_imagefile) as f:
    v_img = f[0].data

vr_transform, (r_locs, v_locs1) = aa.find_transform(r_img, v_img)
vb_transform, (b_locs, v_locs2) = aa.find_transform(b_img, v_img)

compiled_data = {'X': list(), 'Y': list(), 'R Flux': list(), 'B Flux': list(), 'V Flux': list()}
for _, row in v_fluxes.iterrows():
    loc = (row['X'], row['Y'])
    r_locT = aa.matrix_transform(loc, vr_transform.params)
    r_locT = (r_locT[0, 0], r_locT[0, 1])
    seps2 = [distance2(loc, vloc2) for vloc2 in v_locs2]
    b_locT = aa.matrix_transform(loc, vb_transform.params)
    b_locT = (b_locT[0, 0], b_locT[0, 1])
    b_flux = []
    r_flux = []
    for _, b_row in b_fluxes.iterrows():
        b_loc = (b_row['X'], b_row['Y'])
        if distance2(b_loc, b_locT) < FWHM:
            b_flux.append(b_row['Flux'])
    for _, r_row in r_fluxes.iterrows():
        r_loc = (r_row['X'], r_row['Y'])
        if distance2(r_loc, r_locT) < FWHM:
            r_flux.append(r_row['Flux'])
    if len(b_flux) == 1 and len(r_flux) == 1:
        compiled_data['X'].append(loc[0])
        compiled_data['Y'].append(loc[1])
        compiled_data['R Flux'].append(r_flux[0])
        compiled_data['B Flux'].append(b_flux[0])
        compiled_data['V Flux'].append(row['Flux'])

pd.DataFrame(compiled_data).to_csv(output_file, index=False)
