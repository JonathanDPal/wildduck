import numpy as np
import pandas as pd
import sys
from scipy.stats import linregress
from scipy.signal import savgol_filter

# USER INPUTS #
_, ref_direc, target_direc = sys.argv
ref_file, target_file = f'{ref_direc}/{ref_direc}-fluxes.csv', f'{target_direc}/{target_direc}-fluxes.csv'
# END OF USER INPUTS #d

distance_ratios = list()
for filter1, filter2 in zip(['B', 'B', 'V'], ['V', 'R', 'R']):
    ref_fluxes = pd.read_csv(ref_file).sort_values('V Flux', ascending=False)
    median_flux1 = ref_fluxes['V Flux'].median()
    ref_fluxes['Magnitude'] = [-2.5 * np.log10(flux / median_flux1) for flux in ref_fluxes['V Flux']]
    ref_fluxes['Color'] = [-2.5 * np.log10(row[f'{filter1} Flux'] / row[f'{filter2} Flux'])
                           for _, row in ref_fluxes.iterrows()]
    relevant_fluxes1 = ref_fluxes[(-2 < ref_fluxes['Magnitude']) & (ref_fluxes['Magnitude'] < 2)
                                  & (0 < ref_fluxes['Color']) & (ref_fluxes['Color'] < 2)]
    linfit1 = linregress(relevant_fluxes1['Color'], relevant_fluxes1['Magnitude'])
    slope1, intercept1 = linfit1.slope, linfit1.intercept
    ref_magsfit = np.array([slope1 * x + intercept1 for x in np.linspace(0, 1.5, 50000)])
    ref_fluxesfit = 100 ** (-1 * ref_magsfit / 5) * median_flux1

    target_fluxes = pd.read_csv(target_file).sort_values('V Flux', ascending=False)
    median_flux2 = target_fluxes['V Flux'].median()
    target_fluxes['Magnitude'] = [-2.5 * np.log10(flux / median_flux2) for flux in target_fluxes['V Flux']]
    target_fluxes['Color'] = [-2.5 * np.log10(row[f'{filter1} Flux'] / row[f'{filter2} Flux'])
                              for _, row in target_fluxes.iterrows()]
    relevant_fluxes2 = target_fluxes[(-2 < target_fluxes['Magnitude']) & (target_fluxes['Magnitude'] < 2)
                                     & (0 < target_fluxes['Color']) & (target_fluxes['Color'] < 2)]
    linfit2 = linregress(relevant_fluxes2['Color'], relevant_fluxes2['Magnitude'])
    slope2, intercept2 = linfit2.slope, linfit2.intercept
    target_magsfit = np.array([slope2 * x + intercept2 for x in np.linspace(0, 1.5, 50000)])
    target_fluxesfit = 100 ** (-1 * target_magsfit / 5) * median_flux2

    dratios = np.sqrt(ref_fluxesfit / target_fluxesfit)
    distance_ratios.append(np.mean(dratios))

print(distance_ratios)
print(np.mean(distance_ratios))