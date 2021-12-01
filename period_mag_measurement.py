import numpy as np
import pandas as pd
import sys

magfile = 'CY_Aqr_Data/Becca_CY_Aqr_9-28.csv'
df = pd.read_csv(magfile)
df = df[df['Magnitude'] < -0.71]
x, y = df['JD'], df['Magnitude']

cp1 = 1
while np.min(y[:cp1]) > np.min(y[:cp1 + 150]):
    cp1 += 1
cp1 -= 1

cp2 = cp1 + 76
while np.max(y[cp1+75:cp2]) < np.max(y[cp1+75:cp2 + 150]):
    cp2 += 1
cp2 -= 1

cp3 = cp2 + 76
while np.min(y[cp2+75:cp3]) > np.min(y[cp2+75:cp3 + 150]):
    cp3 += 1
cp3 -= 1

cp4 = cp3 + 76
while np.max(y[cp3+75: cp4]) < np.max(y[cp3 + 75: cp4 + 150]):
    cp4 += 1
cp4 -= 1

period = np.mean([x[cp4] - x[cp2], x[cp3] - x[cp1]]) * 24 * 60
magnitude_change = np.mean([y[cp4], y[cp2]]) - np.mean([y[cp3], y[cp1]])
minutes = np.floor(period)
seconds = np.round((period - minutes) * 60)

print(f'Period: {int(minutes)} Minutes, {int(seconds)} Seconds')
print(f'Amplitude: {round(magnitude_change, 3)}')