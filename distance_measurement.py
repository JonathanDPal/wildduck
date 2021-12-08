import numpy as np
import pandas as pd
import sys
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# USER INPUTS #
_, ref_direc, target_direc, fitdeg = sys.argv
ref_file, target_file = f'{ref_direc}/{ref_direc}-fluxes.csv', f'{target_direc}/{target_direc}-fluxes.csv'
# END OF USER INPUTS #


def deg2poly(x, a, b, c):
    return a * x ** 2 + b * x + c


def deg3poly(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d


distance_ratios = list()
for filter1, filter2 in zip(['B', 'B', 'V'], ['V', 'R', 'R']):
    ref_fluxes = pd.read_csv(ref_file).sort_values('V Flux', ascending=False)
    median_flux1 = ref_fluxes['V Flux'].median()
    ref_fluxes['Magnitude'] = [-2.5 * np.log10(flux / median_flux1) for flux in ref_fluxes['V Flux']]
    ref_fluxes['Color'] = [-2.5 * np.log10(row[f'{filter1} Flux'] / row[f'{filter2} Flux'])
                           for _, row in ref_fluxes.iterrows()]
    relevant_fluxes1 = ref_fluxes[(-2 < ref_fluxes['Magnitude']) & (ref_fluxes['Magnitude'] < 2)
                                  & (0 < ref_fluxes['Color']) & (ref_fluxes['Color'] < 2)]

    target_fluxes = pd.read_csv(target_file).sort_values('V Flux', ascending=False)
    median_flux2 = target_fluxes['V Flux'].median()
    target_fluxes['Magnitude'] = [-2.5 * np.log10(flux / median_flux2) for flux in target_fluxes['V Flux']]
    target_fluxes['Color'] = [-2.5 * np.log10(row[f'{filter1} Flux'] / row[f'{filter2} Flux'])
                              for _, row in target_fluxes.iterrows()]
    relevant_fluxes2 = target_fluxes[(-2 < target_fluxes['Magnitude']) & (target_fluxes['Magnitude'] < 2)
                                     & (0 < target_fluxes['Color']) & (target_fluxes['Color'] < 2)]

    allcolors = list(relevant_fluxes1['Color']) + list(relevant_fluxes2['Color'])
    xmin, xmax = np.percentile(allcolors, (10, 90))

    if fitdeg == 'one':
        linfit1 = linregress(relevant_fluxes1['Color'], relevant_fluxes1['Magnitude'])
        slope1, intercept1 = linfit1.slope, linfit1.intercept
        ref_magsfit = np.array([slope1 * x + intercept1 for x in np.linspace(xmin, xmax, 50000)])
        linfit2 = linregress(relevant_fluxes2['Color'], relevant_fluxes2['Magnitude'])
        slope2, intercept2 = linfit2.slope, linfit2.intercept
        target_magsfit = np.array([slope2 * x + intercept2 for x in np.linspace(xmin, xmax, 50000)])
    elif fitdeg == 'two':
        (A1, B1, C1), _ = curve_fit(deg2poly, relevant_fluxes1['Color'], relevant_fluxes1['Magnitude'])
        ref_magsfit = np.array([deg2poly(x, A1, B1, C1) for x in np.linspace(xmin, xmax, 50000)])
        (A2, B2, C2), _ = curve_fit(deg2poly, relevant_fluxes2['Color'], relevant_fluxes2['Magnitude'])
        target_magsfit = np.array([deg2poly(x, A2, B2, C2) for x in np.linspace(xmin, xmax, 50000)])
    elif fitdeg == 'three':
        (A1, B1, C1, D1), _ = curve_fit(deg3poly, relevant_fluxes1['Color'], relevant_fluxes1['Magnitude'])
        ref_magsfit = np.array([deg3poly(x, A1, B1, C1, D1) for x in np.linspace(xmin, xmax, 50000)])
        (A2, B2, C2, D2), _ = curve_fit(deg3poly, relevant_fluxes2['Color'], relevant_fluxes2['Magnitude'])
        target_magsfit = np.array([deg3poly(x, A2, B2, C2, D2) for x in np.linspace(xmin, xmax, 50000)])

    ref_fluxesfit = 100 ** (-1 * ref_magsfit / 5) * median_flux1
    target_fluxesfit = 100 ** (-1 * target_magsfit / 5) * median_flux2
    dratios = np.sqrt(ref_fluxesfit / target_fluxesfit)
    distance_ratios.append(dratios)

print(np.mean(np.array(distance_ratios).flatten()), np.std(np.array(distance_ratios).flatten()))
