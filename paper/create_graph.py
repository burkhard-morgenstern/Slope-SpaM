import csv
import math
import sys

import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(9, 9))

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xticks(np.arange(0, 1.01, 0.1))
plt.yticks(np.arange(0, 1.01, 0.1))
plt.grid(True)

plt.xlabel('distance')
plt.ylabel('estimated distance')

plt.plot([0, 1], [0, 1], color='0000', linewidth=1)

labels = [
    'assembled/assembled',
    'assembled/reads',
    'reads/reads'
]

colors = [
    'C0',
    'C1',
    'C2'
]

for file in sys.argv[1:]:
#for file, label, color in zip(sys.argv[1:], labels, colors):
    rows = []
    with open(file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in reader:
            rows.append(row)
    rows = sorted(rows, key = lambda row: float(row[0]))
    keys = [float(row[0]) for row in rows]
    expected = [float(row[1]) for row in rows]
    values = [list(map(float, row[2:])) for row in rows]
    values = [list(filter(lambda x: not math.isnan(x), l)) for l in values]

    x = np.array(keys)
    y = np.array([np.average(l) for l in values])
    e = np.array([np.std(l) for l in values])

    plt.errorbar(x, y, e, marker='o', capsize=2)
    #plt.errorbar(x, y, e, marker='o', capsize=2, label=label, color=color)

#plt.legend(loc='upper left')
plt.savefig('20000_length_10000.pdf')
plt.show()
