from astropy.io import fits
from matplotlib import pyplot as plt
import sys
import numpy as np

_, fits_file, output_file = sys.argv
with fits.open(fits_file) as f:
    image_array = f[0].data
plt.figure(figsize=(10, 10))
plt.imshow(image_array, vmin = np.percentile(image_array, 40), vmax = np.percentile(image_array, 99.6))
plt.savefig(output_file)
