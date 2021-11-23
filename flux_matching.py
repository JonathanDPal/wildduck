import numpy as np
import pandas as pd

# USER INPUTS #
# CSV files that were output of 'flux_measurement.py'
r_fluxtable = ''
b_fluxtable = ''
v_fluxtable = ''
FWHM = 3
num_sat_stars = 3
output_file = ''
# END OF USER INPUTS #


def distance2(loc1, loc2):
    """
    Distance squared
    """
    x0, y0, _ = loc1
    x1, y1, _ = loc2
    return (x1 - x0) ** 2 + (y1 - y0) ** 2


r_fluxes = pd.read_csv(r_fluxtable).sort_values('Flux', ascending=True)[:num_sat_stars]
b_fluxes = pd.read_csv(b_fluxtable).sort_values('Flux', ascending=True)[:num_sat_stars]
v_fluxes = pd.read_csv(v_fluxtable).sort_values('Flux', ascending=True)[:num_sat_stars]

compiled_data = {'X': list(), 'Y': list(), 'R Flux': list(), 'B Flux': list(), 'V Flux': list()}
for _, row in r_fluxes.iterrows():
    loc = (row['X'], row['Y'])
    b_flux = []
    v_flux = []
    for _, b_row in b_fluxes.iterrows():
        b_loc = (b_row['X'], b_row['Y'])
        if distance2(b_loc, loc) < FWHM:
            b_flux.append(b_row['Flux'])
    for _, v_row in v_fluxes.iterrows():
        v_loc = (v_row['X'], v_row['Y'])
        if distance2(v_loc, loc) < FWHM:
            v_flux.append(v_row['Flux'])
    if len(b_flux) == 1 and len(v_flux) == 1:
        compiled_data['X'].append(loc[0])
        compiled_data['Y'].append(loc[1])
        compiled_data['R Flux'].append(row['Flux'])
        compiled_data['B Flux'].append(b_flux[0])
        compiled_data['V Flux'].append(v_flux[0])

pd.DataFrame(compiled_data).to_csv('matchedfluxes.csv', index=False)
