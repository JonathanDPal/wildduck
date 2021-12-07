import numpy as np
import pandas as pd
import sys
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from matplotlib import pyplot as plt

_, direc, fittype = sys.argv

ref_file = 'm52/m52-fluxes.csv'
filter1 = 'B'
filter2 = 'V'


def deg2poly(x, a, b, c):
    return a * x ** 2 + b * x + c


def deg3poly(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d


ref_fluxes = pd.read_csv(ref_file).sort_values('V Flux', ascending=False)
median_flux1 = ref_fluxes['V Flux'].median()
ref_fluxes['Magnitude'] = [-2.5 * np.log10(flux / median_flux1) for flux in ref_fluxes['V Flux']]
ref_fluxes['Color'] = [-2.5 * np.log10(row[f'{filter1} Flux'] / row[f'{filter2} Flux'])
                       for _, row in ref_fluxes.iterrows()]
relevant_fluxes1 = ref_fluxes[(-2 < ref_fluxes['Magnitude']) & (ref_fluxes['Magnitude'] < 2) &
                              (0 < ref_fluxes['Color']) & (ref_fluxes['Color'] < 2)]
if fittype == 'one':
    linfit1 = linregress(relevant_fluxes1['Color'], relevant_fluxes1['Magnitude'])
    slope1, intercept1 = linfit1.slope, linfit1.intercept
    ref_magsfit = np.array([slope1 * x + intercept1 for x in np.linspace(0, 1.5, 50000)])
elif fittype == 'two':
    (A, B, C), _ = curve_fit(deg2poly, relevant_fluxes1['Color'], relevant_fluxes1['Magnitude'])
    ref_magsfit = np.array([deg2poly(x, A, B, C) for x in np.linspace(0, 1.5, 50000)])
elif fittype == 'three':
    (A, B, C, D), _ = curve_fit(deg3poly, relevant_fluxes1['Color'], relevant_fluxes1['Magnitude'])
    ref_magsfit = np.array([deg3poly(x, A, B, C, D) for x in np.linspace(0, 1.5, 50000)])
relevant_fluxes1.plot.scatter('Color', 'Magnitude', c='red')
plt.plot(np.linspace(0, 1.5, 50000), ref_magsfit, color='blue')
plt.show()
