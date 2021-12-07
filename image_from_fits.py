from astropy.io import fits
from PIL import Image
import sys

_, fits_file, img_output = sys.argv
with fits.open(fits_file) as f:
    image_array = f[0].data
im = Image.fromarray(image_array)
im = im.convert('RGB')
im.save(img_output)
