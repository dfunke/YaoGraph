import os
import re

import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

reFilename = re.compile("benchmark_([a-zA-Z0-9]+)_([a-zA-Z]+).csv")
isRunningInPyCharm = "PYCHARM_HOSTED" in os.environ

DIR = '/home/funke/devel/geograph/benchmark/data'
EXT = 'png'


def saveFig(f: str):
    plt.tight_layout()
    plt.savefig(f)
    if isRunningInPyCharm:
        plt.show()
    plt.close()


sns.set_style('whitegrid')

######################################################################################################
# %% Priority Queue plots
pqDataFile = os.path.join(DIR, 'pq.csv')
if os.path.exists(pqDataFile):
    with open(pqDataFile, 'r') as file:
        header = file.readline()

    header = header[1:].strip().split()

    gPQ = pd.read_csv(pqDataFile, sep=' ', comment='#', names=header)
    lPQ = gPQ[gPQ['k'] == 0]

    plt.stackplot(lPQ['step'], lPQ['ipPro'], lPQ['isPro'], lPQ['delPro'],
                  labels=['input points', 'intersections', 'deletions'])
    plt.legend(loc='upper left')
    plt.xlabel('algorithm step')
    plt.ylabel('events processed')
    saveFig("60_pq_events_processed.%s" % (EXT))

    plt.stackplot(lPQ['step'], lPQ['ipQ'], lPQ['isQ'], lPQ['delQ'],
                  labels=['input points', 'intersections', 'deletions'])
    plt.legend()
    plt.xlabel('algorithm step')
    plt.ylabel('events in PQ')
    saveFig("60_pq_events_inQueue.%s" % (EXT))

    plt.stackplot(lPQ['step'], lPQ['isQ'], lPQ['delQ'], labels=['intersections', 'deletions'])
    plt.legend()
    plt.xlabel('algorithm step')
    plt.ylabel('events in PQ')
    saveFig("60_pq_events_inQueue_noInput.%s" % (EXT))

    plt.stackplot(lPQ['step'], lPQ['slSize'], labels=['rays'])
    plt.legend()
    plt.xlabel('algorithm step')
    plt.ylabel('rays in SL')
    saveFig("60_rays_inSL.%s" % (EXT))

######################################################################################################
# %% Priority Queue aggregate plots
statsDataFile = os.path.join(DIR, 'stats.csv')
if os.path.exists(statsDataFile):
    with open(statsDataFile, 'r') as file:
        header = file.readline()

    header = header[1:].strip().split()

    gStats = pd.read_csv(statsDataFile, sep=' ', comment='#', names=header)
    gStats['maxIsDelQ'] = gStats['maxIsQ'] + gStats['maxDelQ']
    gStats['stepsPN'] = gStats['steps'] / gStats['n']

    sns.lineplot(data=gStats, x='n', y='stepsPN', hue='dist', style='dist', markers=True)
    plt.legend()
    plt.xlabel('n')
    plt.ylabel('events/n')
    saveFig("70_pq_events_processed.%s" % (EXT))

    sns.lineplot(data=gStats, x='n', y='maxIsDelQ', hue='dist', style='dist', markers=True)
    plt.legend()
    plt.xlabel('n')
    plt.ylabel('events in PQ')
    saveFig("70_pq_events_inQueue_noInput.%s" % (EXT))

    sns.lineplot(data=gStats, x='n', y='maxSlSize', hue='dist', style='dist', markers=True)
    plt.legend()
    plt.xlabel('n')
    plt.ylabel('rays in SL')
    saveFig("70_rays_inSL.%s" % (EXT))

######################################################################################################
# %% Sweepline tree update operations plots

upDataFile = os.path.join(DIR, 'update.csv')
if os.path.exists(upDataFile):
    with open(upDataFile, 'r') as file:
        header = file.readline()

    header = header[1:].strip().split()

    gUp = pd.read_csv(upDataFile, sep=' ', comment='#', names=header)
    sns.displot(gUp, x='change', col='op', discrete=True, stat='probability', common_norm=False, shrink=.8)
    plt.xticks([0, 1])

    saveFig("50_tree_updates.%s" % (EXT))

######################################################################################################
# %% Runtime plots - read data

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

gData['tpn'] = gData['t'] / gData['n']

######################################################################################################
# %% Runtime plots - optics

snsPalette = sns.color_palette(n_colors=3)
KernelPalette = {'InexactKernel': snsPalette[0], 'CGALExactPredInexactCon': snsPalette[1],
                 'CGALExactPredExactCon': snsPalette[2]}

snsPalette = sns.color_palette(n_colors=5)
DistPalette = {'gaussian': snsPalette[0], 'uni': snsPalette[1], 'grid': snsPalette[2], 'road': snsPalette[3],
               'star': snsPalette[4]}

snsPalette = sns.color_palette(n_colors=4)
AlgorithmDash = {'NaiveYao': (0, 5, 10), 'Sweepline': (None, None), 'GridYao100': (0, 1, 1), 'CGALYao': (0, 3, 5, 1, 5)}
AlgorithmMarkers = {'NaiveYao': 's', 'Sweepline': 'D', 'GridYao100': 'o', 'CGALYao': '^'}
AlgorithmPalette = {'NaiveYao': snsPalette[0], 'Sweepline': snsPalette[1], 'GridYao100': snsPalette[2],
                    'CGALYao': snsPalette[3]}

figTextX = .94
figTextY = .91
######################################################################################################
# %% Runtime plots - overview plot

sns.lineplot(data=gData, x='n', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash,
             markers=AlgorithmMarkers,
             palette=KernelPalette)
plt.xscale('log')
plt.yscale('log', base=2)
saveFig("01_runtime.%s" % (EXT))

######################################################################################################
# %% Runtime plots - per kernel plots

for kernel in gData['kernel'].unique():
    # filter data by dist
    fData = gData[gData['kernel'] == kernel]

    sns.lineplot(data=fData, x='n', y='tpn', style='algorithm', hue='algorithm', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=AlgorithmPalette,
                 markersize=8)
    plt.xscale('log')
    plt.yscale('log', base=2)
    plt.text(figTextX, figTextY, 'Kernel: %s' % (kernel),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig("01_%s_runtime.%s" % (kernel, EXT))

    sns.lineplot(data=fData, x='n', y='tpn', style='algorithm', hue='dist', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=DistPalette,
                 markersize=8)
    plt.xscale('log')
    plt.yscale('log', base=2)
    plt.text(figTextX, figTextY, 'Kernel: %s' % (kernel),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig("01_%s_runtime_dist.%s" % (kernel, EXT))

######################################################################################################
# %% Runtime plots - per dist plots

for dist in gData['dist'].unique():

    # filter data by dist
    fData = gData[gData['dist'] == dist]

    sns.lineplot(data=fData, x='n', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=KernelPalette)
    plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig("01_%s_runtime.%s" % (dist, EXT))

    for alg in fData['algorithm'].unique():
        sns.lineplot(data=fData[fData['algorithm'] == alg], x='n', y='tpn', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
                 horizontalalignment='right',
                 verticalalignment='bottom',
                 transform=plt.gcf().transFigure)
        saveFig("02_%s_runtime_alg_%s.%s" % (dist, alg, EXT))

    for kernel in fData['kernel'].unique():
        sns.lineplot(data=fData[fData['kernel'] == kernel], x='n', y='tpn', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
                 horizontalalignment='right',
                 verticalalignment='bottom',
                 transform=plt.gcf().transFigure)
        saveFig("03_%s_runtime_kernel_%s.%s" % (dist, kernel, EXT))

    sns.lineplot(data=fData[~fData['algorithm'].isin(['NaiveYao', 'CGALYao'])], x='n', y='tpn', style='algorithm',
                 hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers, palette=KernelPalette)
    plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig("04_%s_runtime_sweepline_grid.%s" % (dist, EXT))

    sns.lineplot(
        data=fData[
            ~fData['algorithm'].isin(['NaiveYao', 'CGALYao']) & ~fData['kernel'].isin(['CGALExactPredExactCon'])],
        x='n', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers,
        palette=KernelPalette)
    plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig("05_%s_runtime_sweepline_grid_noExactCon.%s" % (dist, EXT))

