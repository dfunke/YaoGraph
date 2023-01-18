import os
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pylab as plt
import seaborn as sns

reRuntimeFilename = re.compile('benchmark_([a-zA-Z0-9]+)_([a-zA-Z]+).csv')
reConesFilename = re.compile('benchmark_cones_([a-zA-Z0-9]+)_([a-zA-Z]+).csv')
isRunningInPyCharm = 'PYCHARM_HOSTED' in os.environ

BASE_DIR = '/home/funke/devel/geograph/benchmark/'
DATA_DIR = os.path.join(BASE_DIR, 'data')
PLOT_DIR = os.path.join(BASE_DIR, 'plots')
EXT = 'png'


def saveFig(f: str):
    plt.tight_layout(pad=0.3)
    plt.savefig(os.path.join(PLOT_DIR, f))
    if isRunningInPyCharm:
        plt.show()
    plt.close()

BASE_SIZE = 20
LARGE_SIZE = 34

# sns.set(style='whitegrid', font_scale=1, rc={'text.usetex': True,
#                                                'text.latex.preamble': r'\usepackage{nicefrac} \usepackage{amsmath} \usepackage{lmodern}',
#                                                "font.size": BASE_SIZE,
#                                              "font.family": 'serif'
#                                              "axes.titlesize": BASE_SIZE,
#                                              "axes.labelsize": LARGE_SIZE,
#                                              "xtick.labelsize" : LARGE_SIZE, "ytick.labelsize" : LARGE_SIZE,
#                                              "legend.fontsize" : BASE_SIZE,
#                                              "legend.title_fontsize" : BASE_SIZE})

sns.set(style='whitegrid', font_scale=2, rc={'text.usetex': True,
                                               'text.latex.preamble': r'\usepackage{nicefrac} \usepackage{amsmath} \usepackage{lmodern}',
                                             "font.family": 'serif'})

LegendLabels = {
    'gaussian': r'Gaussian', 'uni': r'Uniform', 'grid': r'Grid', 'road': r'Road', 'stars': r'Stars',
    'NaiveYao': r'Naive Yao', 'Sweepline': r'Sweepline', 'GridYao': r'Grid Yao', 'CGALYao': r'CGAL Yao',
    'InexactKernel': r'Inexact Kernel', 'CGALExactPredInexactCon': r'CGAL EPIC', 'CGALExactPredExactCon': r'CGAL EPEC',
    'kernel': 'Kernel', 'algorithm': r'Algorithm', 'dist': r'Distribution',
}


def setLegend(ax: mpl.axes.Axes):
    L = ax.get_legend()

    if (L.get_title() and L.get_title().get_text() != ''):
        L.set_title(LegendLabels[L.get_title().get_text()])
        L.set_title(None)

    for i in range(len(L.get_texts())):
        L.get_texts()[i].set_text(LegendLabels[L.get_texts()[i].get_text()])


def unsetLegend(ax: mpl.axes.Axes):
    ax.get_legend().remove()


######################################################################################################
# %% Priority Queue plots
pltGroup = 10
pqDataFile = os.path.join(DATA_DIR, 'pq.csv')
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
    saveFig('%i_pq_events_processed.%s' % (pltGroup, EXT))

    plt.stackplot(lPQ['step'], lPQ['ipQ'], lPQ['isQ'], lPQ['delQ'],
                  labels=['input points', 'intersections', 'deletions'])
    # setLegend()
    plt.xlabel('algorithm step')
    plt.ylabel('events in PQ')
    saveFig('%i_pq_events_inQueue.%s' % (pltGroup, EXT))

    plt.stackplot(lPQ['step'], lPQ['isQ'], lPQ['delQ'], labels=['intersections', 'deletions'])
    # setLegend()
    plt.xlabel('algorithm step')
    plt.ylabel('events in PQ')
    saveFig('%i_pq_events_inQueue_noInput.%s' % (pltGroup, EXT))

    plt.stackplot(lPQ['step'], lPQ['slSize'], labels=['rays'])
    # setLegend()
    plt.xlabel('algorithm step')
    plt.ylabel('rays in SL')
    saveFig('%i_rays_inSL.%s' % (pltGroup, EXT))

######################################################################################################
# %% Priority Queue aggregate plots
pltGroup = 20
statsDataFile = os.path.join(DATA_DIR, 'stats.csv')
if os.path.exists(statsDataFile):
    with open(statsDataFile, 'r') as file:
        header = file.readline()

    header = header[1:].strip().split()

    gStats = pd.read_csv(statsDataFile, sep=' ', comment='#', names=header)
    gStats['maxIsDelQ'] = gStats['maxIsQ'] + gStats['maxDelQ']
    gStats['stepsPN'] = gStats['steps'] / gStats['n']
    x = gStats['n'].unique()
    x.sort()
    sqrtN = np.sqrt(x)
    gStats['maxIsDelQSqrt'] = gStats['maxIsDelQ'] / np.sqrt(gStats['n'])
    gStats['maxSlSizeSqrt'] = gStats['maxSlSize'] / np.sqrt(gStats['n'])

    ax = sns.lineplot(data=gStats, x='n', y='stepsPN', hue='dist', style='dist', markers=True)
    unsetLegend(ax)
    plt.xscale('log')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\nicefrac{\text{events}}{n}$')
    saveFig('%i_pq_events_processed.%s' % (pltGroup, EXT))

    ax = sns.lineplot(data=gStats, x='n', y='maxIsDelQSqrt', hue='dist', style='dist', markers=True)
    unsetLegend(ax)
    plt.xscale('log')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\nicefrac{\text{events in PQ}}{\sqrt{n}}$')
    saveFig('%i_pq_events_inQueue_noInput.%s' % (pltGroup, EXT))

    ax = sns.lineplot(data=gStats, x='n', y='maxSlSizeSqrt', hue='dist', style='dist', markers=True)
    setLegend(ax)
    plt.xscale('log')
    plt.xlabel('$n$')
    plt.ylabel(r'$\nicefrac{\text{rays in SL}}{\sqrt{n}}$')
    saveFig('%i_rays_inSL.%s' % (pltGroup, EXT))

######################################################################################################
# %% Sweepline tree update operations plots
pltGroup = 30
upDataFile = os.path.join(DATA_DIR, 'update.csv')
if os.path.exists(upDataFile):
    with open(upDataFile, 'r') as file:
        header = file.readline()

    header = header[1:].strip().split()

    gUp = pd.read_csv(upDataFile, sep=' ', comment='#', names=header)
    sns.displot(gUp, x='change', col='op', discrete=True, stat='probability', common_norm=False, shrink=.8)
    plt.xticks([0, 1])

    saveFig('%i_tree_updates.%s' % (pltGroup, EXT))

######################################################################################################
# %% Runtime plots - read data

llData = []
for f in os.listdir(DATA_DIR):
    match = reRuntimeFilename.match(f)
    if match:
        # get header from first commented line
        with open(os.path.join(DATA_DIR, f), 'r') as file:
            header = file.readline()

        header = header[1:].strip().split()

        lData = pd.read_csv(os.path.join(DATA_DIR, f), sep=' ', comment='#', names=header)
        lData['algorithm'] = match.group(1)
        lData['kernel'] = match.group(2)
        llData.append(lData)

gData = pd.concat(llData)
gData.set_index(['algorithm', 'kernel', 'dist', 'n', 'seed', 'rep'], inplace=True)
gData.reset_index(inplace=True)

gData['tpn'] = gData['t'] / gData['n']
gData['speedup'] = 1

# calculate speedup over naive
i = 1
for dist in gData['dist'].unique():
    for alg in gData['algorithm'].unique():
        for kernel in gData[gData['algorithm'] == alg]['kernel'].unique():
            maskAlg = (gData['dist'] == dist) & (gData['algorithm'] == alg) & (gData['kernel'] == kernel)
            maskNav = (gData['dist'] == dist) & (gData['algorithm'] == 'NaiveYao') & (gData['kernel'] == kernel)

            gData.loc[maskAlg, 'speedup'] = gData.loc[maskNav, 't'].array / gData.loc[maskAlg, 't']

######################################################################################################
# %% Runtime plots - optics

snsPalette = sns.color_palette(n_colors=3)
KernelPalette = {'InexactKernel': snsPalette[0], 'CGALExactPredInexactCon': snsPalette[1],
                 'CGALExactPredExactCon': snsPalette[2]}
KernelDash = {'InexactKernel': (0, 5, 10), 'CGALExactPredInexactCon': (None, None), 'CGALExactPredExactCon': (0, 1, 1)}
KernelMarkers = {'InexactKernel': 's', 'CGALExactPredInexactCon': 'D', 'CGALExactPredExactCon': 'o'}

snsPalette = sns.color_palette(n_colors=5)
DistPalette = {'gaussian': snsPalette[0], 'uni': snsPalette[1], 'grid': snsPalette[2], 'road': snsPalette[3],
               'stars': snsPalette[4]}

snsPalette = sns.color_palette(n_colors=4)
AlgorithmDash = {'NaiveYao': (0, 5, 10), 'Sweepline': (None, None), 'GridYao': (0, 1, 1), 'CGALYao': (0, 3, 5, 1, 5)}
AlgorithmMarkers = {'NaiveYao': 's', 'Sweepline': 'D', 'GridYao': 'o', 'CGALYao': '^'}
AlgorithmPalette = {'NaiveYao': snsPalette[0], 'Sweepline': snsPalette[1], 'GridYao': snsPalette[2],
                    'CGALYao': snsPalette[3]}

figTextX = .94
figTextY = .91

######################################################################################################
# %% Runtime plots - overview plot

pltGroup = 40
ax = sns.lineplot(data=gData, x='n', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash,
             markers=AlgorithmMarkers,
             palette=KernelPalette)
w, h = plt.gcf().get_size_inches()
plt.gcf().set_size_inches(1.5*w, h)
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.1))
setLegend(ax)
plt.xscale('log')
plt.yscale('log', base=2)
plt.xlabel(r'$n$')
plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
saveFig('%i_runtime.%s' % (pltGroup, EXT))

ax = sns.lineplot(data=gData, x='n', y='speedup', style='kernel', hue='algorithm', dashes=KernelDash,
             markers=KernelMarkers,
             palette=AlgorithmPalette)
w, h = plt.gcf().get_size_inches()
plt.gcf().set_size_inches(1.5*w, h)
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.1))
setLegend(ax)
plt.xscale('log')
plt.yscale('log', base=2)
plt.xlabel(r'$n$')
plt.ylabel(r'speedup')
saveFig('%i_speedup.%s' % (pltGroup, EXT))

######################################################################################################
# %% Runtime plots - overview plot (no Naive/CGAL)
pltGroup = 45
fData = gData[(gData['algorithm'] == 'Sweepline') ^ (gData['algorithm'] == 'GridYao')]
ax = sns.lineplot(data=fData, x='n', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash,
             markers=AlgorithmMarkers,
             palette=KernelPalette)
unsetLegend(ax)
plt.xscale('log')
plt.yscale('log', base=2)
plt.xlabel(r'$n$')
plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
saveFig('%i_runtime.%s' % (pltGroup, EXT))

ax = sns.lineplot(data=fData, x='n', y='speedup', style='kernel', hue='algorithm', dashes=KernelDash,
             markers=KernelMarkers,
             palette=AlgorithmPalette)
unsetLegend(ax)
plt.xscale('log')
plt.yscale('log', base=2)
plt.xlabel(r'$n$')
plt.ylabel(r'speedup')
saveFig('%i_speedup.%s' % (pltGroup, EXT))

######################################################################################################
# %% Runtime plots - per kernel plots
pltGroup = 50
for kernel in gData['kernel'].unique():
    # filter data by dist
    fData = gData[gData['kernel'] == kernel]

    ax = sns.lineplot(data=fData, x='n', y='tpn', style='algorithm', hue='algorithm', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=AlgorithmPalette,
                 markersize=8)
    setLegend(ax)
    plt.xscale('log')
    plt.yscale('log', base=2)
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
    # plt.text(figTextX, figTextY, 'Kernel: %s' % (LegendLabels[kernel]),
    #          horizontalalignment='right',
    #          verticalalignment='bottom',
    #          transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime.%s' % (pltGroup, kernel, EXT))

    sns.lineplot(data=fData, x='n', y='tpn', style='algorithm', hue='dist', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=DistPalette,
                 markersize=8)
    plt.xscale('log')
    plt.yscale('log', base=2)
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
    plt.text(figTextX, figTextY, 'Kernel: %s' % (LegendLabels[kernel]),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime_dist.%s' % (pltGroup, kernel, EXT))

######################################################################################################
# %% Runtime plots - per kernel plots (no Naive/CGAL)
pltGroup = 55

fData = gData[(gData['algorithm'] == 'Sweepline') ^ (gData['algorithm'] == 'GridYao')]
for kernel in gData['kernel'].unique():
    # filter data by dist
    ffData = fData[fData['kernel'] == kernel]

    sns.lineplot(data=ffData, x='n', y='tpn', style='algorithm', hue='algorithm', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=AlgorithmPalette,
                 markersize=8)
    plt.xscale('log')
    plt.yscale('log', base=2)
    plt.text(figTextX, figTextY, 'Kernel: %s' % (kernel),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime.%s' % (pltGroup, kernel, EXT))

    sns.lineplot(data=ffData, x='n', y='tpn', style='algorithm', hue='dist', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=DistPalette,
                 markersize=8)
    plt.xscale('log')
    plt.yscale('log', base=2)
    plt.text(figTextX, figTextY, 'Kernel: %s' % (kernel),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime_dist.%s' % (pltGroup, kernel, EXT))

######################################################################################################
# %% Runtime plots - per dist plots
pltGroup = 60
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
    saveFig('%i_%s_runtime.%s' % (pltGroup, dist, EXT))

    for alg in fData['algorithm'].unique():
        sns.lineplot(data=fData[fData['algorithm'] == alg], x='n', y='tpn', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
                 horizontalalignment='right',
                 verticalalignment='bottom',
                 transform=plt.gcf().transFigure)
        saveFig('%i_%s_runtime_alg_%s.%s' % (pltGroup, dist, alg, EXT))

    for kernel in fData['kernel'].unique():
        sns.lineplot(data=fData[fData['kernel'] == kernel], x='n', y='tpn', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
                 horizontalalignment='right',
                 verticalalignment='bottom',
                 transform=plt.gcf().transFigure)
        saveFig('%i_%s_runtime_kernel_%s.%s' % (pltGroup, dist, kernel, EXT))

    sns.lineplot(data=fData[~fData['algorithm'].isin(['NaiveYao', 'CGALYao'])], x='n', y='tpn', style='algorithm',
                 hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers, palette=KernelPalette)
    plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime_sweepline_grid.%s' % (pltGroup, dist, EXT))

    sns.lineplot(
        data=fData[
            ~fData['algorithm'].isin(['NaiveYao', 'CGALYao']) & ~fData['kernel'].isin(['CGALExactPredExactCon'])],
        x='n', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers,
        palette=KernelPalette)
    plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime_sweepline_grid_noExactCon.%s' % (pltGroup, dist, EXT))

######################################################################################################
# %% Runtime over cones plots - read data

llData = []
for f in os.listdir(DATA_DIR):
    match = reConesFilename.match(f)
    if match:
        # get header from first commented line
        with open(os.path.join(DATA_DIR, f), 'r') as file:
            header = file.readline()

        header = header[1:].strip().split()

        lData = pd.read_csv(os.path.join(DATA_DIR, f), sep=' ', comment='#', names=header)
        lData['algorithm'] = match.group(1)
        lData['kernel'] = match.group(2)
        llData.append(lData)

gData = pd.concat(llData)
gData.set_index(['algorithm', 'kernel', 'dist', 'n', 'seed', 'rep'], inplace=True)
gData.reset_index(inplace=True)

gData['tpn'] = gData['t'] / gData['n']
gData['tpk'] = gData['t'] / gData['k']

######################################################################################################
# %% Runtime over cones plots - overview plot

pltGroup = 70
ax = sns.lineplot(data=gData, x='k', y='tpk', style='algorithm', hue='kernel', dashes=AlgorithmDash,
             markers=AlgorithmMarkers,
             palette=KernelPalette)
plt.legend(loc='upper left', bbox_to_anchor=(.8, 1.1))
setLegend(ax)
# plt.xscale('log')
plt.yscale('log', base=10)
plt.xlabel(r'$k$')
plt.ylabel(r'$\nicefrac{t}{k}$ [ms]')
saveFig('%i_runtime_cones.%s' % (pltGroup, EXT))

