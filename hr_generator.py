import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# USER INPUTS #
input_file = 'wild_duck/wild_duck-fluxes.csv'
filter1 = 'B'
filter2 = 'V'
output_image = 'wild_duck/wild_duck-HRdiagram.png'
# END OF USER INPUTS #

fluxes = pd.read_csv(input_file).sort_values('V Flux', ascending=False)
fluxes = fluxes[:int(np.ceil(len(fluxes) / 5))]
ref_flux = fluxes['V Flux'].median()
fluxes['Magnitude'] = [-2.5 * np.log10(flux / ref_flux) for flux in fluxes['V Flux']]
fluxes['Color'] = [-2.5 * np.log10(row[f'{filter1} Flux'] / row[f'{filter2} Flux']) for _, row in fluxes.iterrows()]
fluxes.plot(x='Color', y='Magnitude', kind='scatter', use_index=False, ylim=(-1, 0.5))
plt.xlabel(f'{filter1} - {filter2} (mag)')
plt.ylabel('Relative Visual Magnitude')
plt.savefig(output_image)
