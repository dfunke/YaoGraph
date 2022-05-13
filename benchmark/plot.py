import os
import re

import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

reFilename = re.compile("benchmark_([a-zA-Z0-9]+)_([a-zA-Z]+).csv")

DIR = '/Users/funke/Development/GeoGraph/benchmark/data'
EXT = 'png'

KernelPalette = {'InexactKernel': 'r', 'CGALExactPredInexactCon': 'g', 'CGALExactPredExactCon': 'b'}
AlgorithmDash = {'NaiveYao': (0, 1, 1), 'Sweepline': (None, None), 'GridYao100': (0, 5, 10), 'CGALYao':   (0, 3, 5, 1, 5)}
AlgorithmMarkers = {'NaiveYao': 'o', 'Sweepline': 'D', 'GridYao100': '*', 'CGALYao': '^'}

llData = []
for f in os.listdir(DIR):
    match = reFilename.match(f)
    if match:
        lData = pd.read_csv(os.path.join(DIR, f), sep=' ')
        lData['algorithm'] = match.group(1)
        lData['kernel'] = match.group(2)
        llData.append(lData)

gData = pd.concat(llData)
# gData.set_index(['algorithm', 'kernel', 'dist', 'n', 'seed', 'rep'], inplace=True)

sns.lineplot(data=gData, x='n', y='t', style='algorithm', hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers,
             palette=KernelPalette)
plt.savefig("01_runtime_all.%s" % EXT)
plt.close()

for alg in gData['algorithm'].unique():
    sns.lineplot(data=gData[gData['algorithm'] == alg], x='n', y='t', style='algorithm', hue='kernel',
                 dashes=AlgorithmDash,
                 markers=AlgorithmMarkers, palette=KernelPalette)
    plt.savefig("02_runtime_alg_%s.%s" % (alg, EXT))
    plt.close()

for kernel in gData['kernel'].unique():
    sns.lineplot(data=gData[gData['kernel'] == kernel], x='n', y='t', style='algorithm', hue='kernel',
                 dashes=AlgorithmDash,
                 markers=AlgorithmMarkers, palette=KernelPalette)
    plt.savefig("03_runtime_kernel_%s.%s" % (kernel, EXT))
    plt.close()

sns.lineplot(data=gData[~gData['algorithm'].isin(['NaiveYao', 'CGALYao'])], x='n', y='t', style='algorithm',
             hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers, palette=KernelPalette)
plt.savefig("04_runtime_sweepline_grid.png")
plt.close()

sns.lineplot(
    data=gData[~gData['algorithm'].isin(['NaiveYao', 'CGALYao']) & ~gData['kernel'].isin(['CGALExactPredExactCon'])],
    x='n', y='t', style='algorithm', hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers,
    palette=KernelPalette)
plt.savefig("05_runtime_sweepline_grid_noExactCon.png")
plt.close()
