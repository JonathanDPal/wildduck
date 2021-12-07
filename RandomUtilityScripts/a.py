import numpy as np
import csv

x, y, z, w, v, r = list(), list(), list(), list(), list(), list()
with open('cyaqrvals.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        x.append(float(row['Timestamp (JD)']))
        y.append(float(row[' Obj1 : Magnitude (Centroid)']))
        z.append(float(row[' Ref1 : Magnitude (Centroid)']))
        w.append(float(row[' Ref2 : Magnitude (Centroid)']))
        v.append(float(row[' Ref3 : Magnitude (Centroid)']))
        r.append(float(row[' Chk1 : Magnitude (Centroid)']))


with open('tolatex.txt', 'w') as f:
    for x0, y0, z0, w0, v0, r0 in zip(x, y, z, w, v, r):
        f.write(f'{x0} & {y0} & {z0} & {w0} & {v0} & {r0} \\\\ \n')
