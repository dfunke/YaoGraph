import os
import re

import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

reFilename = re.compile("benchmark_([a-zA-Z0-9]+)_([a-zA-Z]+).csv")
isRunningInPyCharm = "PYCHARM_HOSTED" in os.environ

DIR = '/home/funke/devel/geograph/benchmark/data'
EXT = 'png'

def saveFig(f : str):
    plt.savefig(f)
    if isRunningInPyCharm:
        plt.show()
    plt.close()

KernelPalette = {'InexactKernel': 'r', 'CGALExactPredInexactCon': 'g', 'CGALExactPredExactCon': 'b'}
AlgorithmDash = {'NaiveYao': (0, 1, 1), 'Sweepline': (None, None), 'GridYao100': (0, 5, 10), 'CGALYao': (0, 3, 5, 1, 5)}
AlgorithmMarkers = {'NaiveYao': 'o', 'Sweepline': 'D', 'GridYao100': '*', 'CGALYao': '^'}


######################################################################################################
# %% Priority Queue plots
pqDataFile = os.path.join(DIR, 'pq.csv')
if os.path.exists(pqDataFile):

    with open(pqDataFile, 'r') as file:
        header = file.readline()

    header = header[1:].strip().split()

    gPQ = pd.read_csv(pqDataFile, sep=' ', comment='#', names=header)
    lPQ = gPQ[gPQ['k'] == 0]

    plt.stackplot(lPQ['step'], lPQ['ipPro'], lPQ['isPro'], lPQ['delPro'], labels=['input points', 'intersections', 'deletions'])
    plt.legend(loc='upper left')
    saveFig("60_pq_events_processed.%s" % (EXT))

    plt.stackplot(lPQ['step'], lPQ['ipQ'], lPQ['isQ'], lPQ['delQ'], labels=['input points', 'intersections', 'deletions'])
    plt.legend()
    saveFig("60_pq_events_inQueue.%s" % (EXT))

    plt.stackplot(lPQ['step'], lPQ['isQ'], lPQ['delQ'], labels=['intersections', 'deletions'])
    plt.legend()
    saveFig("60_pq_events_inQueue_noInput.%s" % (EXT))

    plt.stackplot(lPQ['step'], lPQ['slSize'], labels=['rays'])
    plt.legend()
    saveFig("60_rays_inSL.%s" % (EXT))

######################################################################################################
# %% Sweepline tree update operations plots

upDataFile = os.path.join(DIR, 'update.csv')
if os.path.exists(upDataFile):
    with open(upDataFile, 'r') as file:
        header = file.readline()

    header = header[1:].strip().split()

    gUp = pd.read_csv(upDataFile, sep=' ', comment='#', names=header)
    sns.displot(gUp, x='change', col='op', discrete=True, stat='probability', common_norm=False)

    saveFig("50_tree_updates.%s" % (EXT))

######################################################################################################
# %% Runtime plots

llData = []
for f in os.listdir(DIR):
    match = reFilename.match(f)
    if match:
        # get header from first commented line
        with open(os.path.join(DIR, f), 'r') as file:
            header = file.readline()

        header = header[1:].strip().split()

        lData = pd.read_csv(os.path.join(DIR, f), sep=' ', comment='#', names=header)
        lData['algorithm'] = match.group(1)
        lData['kernel'] = match.group(2)
        llData.append(lData)

gData = pd.concat(llData)
gData.set_index(['algorithm', 'kernel', 'dist', 'n', 'seed', 'rep'], inplace=True)
gData.reset_index(inplace=True)

for dist in gData['dist'].unique():

    # filter data by dist
    fData = gData[gData['dist'] == dist]

    sns.lineplot(data=fData, x='n', y='t', style='algorithm', hue='kernel', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=KernelPalette)
    plt.figtext(.7, .89, 'Distribution: %s' % (dist))
    saveFig("01_%s_runtime.%s" % (dist, EXT))

    for alg in fData['algorithm'].unique():
        sns.lineplot(data=fData[fData['algorithm'] == alg], x='n', y='t', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.figtext(.7, .89, 'Distribution: %s' % (dist))
        saveFig("02_%s_runtime_alg_%s.%s" % (dist, alg, EXT))

    for kernel in fData['kernel'].unique():
        sns.lineplot(data=fData[fData['kernel'] == kernel], x='n', y='t', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.figtext(.7, .89, 'Distribution: %s' % (dist))
        saveFig("03_%s_runtime_kernel_%s.%s" % (dist, kernel, EXT))

    sns.lineplot(data=fData[~fData['algorithm'].isin(['NaiveYao', 'CGALYao'])], x='n', y='t', style='algorithm',
                 hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers, palette=KernelPalette)
    plt.figtext(.7, .89, 'Distribution: %s' % (dist))
    saveFig("04_%s_runtime_sweepline_grid.%s" % (dist, EXT))

    sns.lineplot(
        data=fData[
            ~fData['algorithm'].isin(['NaiveYao', 'CGALYao']) & ~fData['kernel'].isin(['CGALExactPredExactCon'])],
        x='n', y='t', style='algorithm', hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers,
        palette=KernelPalette)
    plt.figtext(.7, .89, 'Distribution: %s' % (dist))
    saveFig("05_%s_runtime_sweepline_grid_noExactCon.%s" % (dist, EXT))
