import numpy as np
from scipy.stats import median_abs_deviation
from glob import glob
from astropy.io import fits
import os
import sys

_, direc = sys.argv
os.chdir(direc)
flatimages = [fits.open(flatfile)[0].data.flatten() for flatfile in glob('flats*')]
darkimages = [fits.open(darkfile)[0].data.flatten() for darkfile in glob('dark*')]
biasimages = [fits.open(biasfile)[0].data.flatten() for biasfile in glob('bias*')]

assert len(np.unique([len(flatimage) for flatimage in flatimages])) == 1
assert len(np.unique([len(darkimage) for darkimage in darkimages])) == 1
assert len(np.unique([len(biasimage) for biasimage in biasimages])) == 1

lfi = len(flatimages[0])
ldi = len(darkimages[0])
lbi = len(biasimages[0])

flaterror = np.mean([median_abs_deviation([flatimage[index] for flatimage in flatimages]) /
                     np.median([flatimage[index] for flatimage in flatimages]) for index in range(lfi)])
darkerror = np.mean([np.std([darkimage[index] for darkimage in darkimages]) /
                     np.mean([darkimage[index] for darkimage in darkimages]) for index in range(ldi)])
biaserror = np.mean([np.std([biasimage[index] for biasimage in biasimages]) /
                     np.mean([biasimage[index] for biasimage in biasimages]) for index in range(lbi)])

with open('errors.txt', 'w') as f:
    f.write(f'{flaterror}\n')
    f.write(f'{darkerror}\n')
    f.write(f'{biaserror}\n')
