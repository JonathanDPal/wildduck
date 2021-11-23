import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# USER INPUTS #
input_file = ''
filter1 = 'B'
filter2 = 'V'
output_image = ''
# END OF USER INPUTS #

fluxes = pd.read_csv(input_file)
ref_flux = fluxes['V Flux'].median()
fluxes['Magnitude'] = [-2.5 * np.log10(flux / ref_flux) for flux in fluxes['V Flux']]
fluxes['Color'] = [-2.5 * np.log10(row['B Flux'] / row['V Flux']) for _, row in fluxes.iterrows()]
fluxes.plot(x='Color', y='Magnitude', kind='scatter', use_index=False)
plt.xlabel('B - V (mag)')
plt.ylabel('Relative Visual Magnitude')
plt.savefig(output_image)
