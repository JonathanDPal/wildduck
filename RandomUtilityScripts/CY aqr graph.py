import numpy as np
import csv
import matplotlib.pyplot as plt

x = []
y = []

with open('cyaqrvals.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        x.append(float(row['Timestamp (JD)']))
        y.append(float(row[' Obj1 : Magnitude (Centroid)']))

y = [100 ** (-y0 / 5) for y0 in y]

plt.scatter(x, y, marker='.', c='red')
plt.scatter((x[18], x[73]), (y[18], y[73]), marker='.', c='blue')
plt.scatter((x[33]), (y[33]), marker='.', c='black')
plt.xlabel('Time (JD)')
plt.ylabel('Relative Flux (ratio)')
plt.title('Light Curve of CY Aqr (V Band)')
plt.yscale('log')
ticks = [4 + 0.5 * n for n in range(6)]
plt.yticks(ticks, [str(tick) for tick in ticks])
plt.savefig('Light Curve.png', dpi=500)

print(24 * (x[73] - x[18]))  # convert JD to hours
print(y[33] - y[18])  # difference in max and min mag
print(y[33] - y[73])
