import numpy as np
import csv

x, y, z, w, v = list(), list(), list(), list(), list()
with open('../WildDuck/0.5FWHM/WildDuck-fluxes.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        x.append(float(row['X']))
        y.append(float(row['Y']))
        z.append(float(row['R Flux']))
        w.append(float(row['B Flux']))
        v.append(float(row['V Flux']))


with open('../tolatex.txt', 'w') as f:
    for x0, y0, z0, w0, v0 in zip(x, y, z, w, v):
        f.write(f'{round(x0, 2)} & {round(y0, 2)} & {round(z0, 2)} & {round(w0, 2)} & {round(v0, 2)} \\\\ \n')
