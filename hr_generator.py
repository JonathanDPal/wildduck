import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys

# USER INPUTS #
_, direc, filter1, filter2 = sys.argv
input_file = f'{direc}/{direc}-fluxes.csv'
output_image = f'{direc}/{direc}-HRdiagram_{filter1}-{filter2}.png'
# END OF USER INPUTS #

fluxes = pd.read_csv(input_file).sort_values('V Flux', ascending=False)
ref_flux = fluxes['V Flux'].median()
fluxes['Magnitude'] = [-2.5 * np.log10(flux / ref_flux) for flux in fluxes['V Flux']]
fluxes['Color'] = [-2.5 * np.log10(row[f'{filter1} Flux'] / row[f'{filter2} Flux']) for _, row in fluxes.iterrows()]
fluxes.plot(x='Color', y='Magnitude', kind='scatter', use_index=False, xlim=(0, 2), ylim=(-2, 2))
plt.xlabel(f'{filter1} - {filter2} (mag)')
plt.ylabel('Relative Visual Magnitude')
plt.savefig(output_image)
