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
                                             "font.family": 'serif',
                                             "lines.markersize": 10})

LegendLabels = {
    'gaussian': r'Gaussian', 'uni': r'Uniform', 'grid': r'Grid', 'road': r'Road', 'stars': r'Stars',
    'circle': r'Circle', 'bubbles': r'Bubbles',
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
# %% Runtime plots - optics

snsPalette = sns.color_palette(n_colors=3)
KernelPalette = {'InexactKernel': snsPalette[0], 'CGALExactPredInexactCon': snsPalette[1],
                 'CGALExactPredExactCon': snsPalette[2]}
KernelDash = {'InexactKernel': (0, 5, 10), 'CGALExactPredInexactCon': (None, None), 'CGALExactPredExactCon': (0, 1, 1)}
KernelMarkers = {'InexactKernel': 's', 'CGALExactPredInexactCon': 'D', 'CGALExactPredExactCon': 'o'}

snsPalette = sns.color_palette(n_colors=7)
DistPalette = {'gaussian': snsPalette[0], 'uni': snsPalette[1], 'grid': snsPalette[2], 'road': snsPalette[3],
               'stars': snsPalette[4], 'circle': snsPalette[5], 'bubbles': snsPalette[6]}

snsPalette = sns.color_palette(n_colors=4)
AlgorithmDash = {'NaiveYao': (0, 5, 10), 'Sweepline': (None, None), 'GridYao': (0, 1, 1), 'CGALYao': (0, 3, 5, 1, 5)}
AlgorithmMarkers = {'NaiveYao': 's', 'Sweepline': 'D', 'GridYao': 'o', 'CGALYao': '^'}
AlgorithmPalette = {'NaiveYao': snsPalette[0], 'Sweepline': snsPalette[1], 'GridYao': snsPalette[2],
                    'CGALYao': snsPalette[3]}

figTextX = .94
figTextY = .91

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
    gStats['maxIsDelQSqrt'] = gStats['maxIsDelQ'] / np.sqrt(gStats['n'])
    gStats['maxIsDelQPN'] = gStats['maxIsDelQ'] / gStats['n']
    gStats['maxSlSizeSqrt'] = gStats['maxSlSize'] / np.sqrt(gStats['n'])
    gStats['maxSlSizePN'] = gStats['maxSlSize'] / (gStats['n'] * np.log(gStats['n']))

    Ns = pd.Series(gStats[gStats['dist'] == 'uni']['n'].unique())


    def custom_round(x):
        return Ns[Ns.sub(x).abs().idxmin()]


    gStats['dN'] = gStats['n'].apply(lambda x: custom_round(x))

    # lStats = gStats[(gStats['dist'] != 'circle') & (gStats['dist'] != 'road') & (gStats['dist'] != 'stars')]
    lStats = gStats[gStats['dist'] != 'bubbles']

    ax = sns.lineplot(data=lStats, x='dN', y='stepsPN', hue='dist', style='dist', markers=True,
                 palette=DistPalette)
    # plt.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0.8))
    unsetLegend(ax)
    plt.xscale('log')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\nicefrac{\text{events}}{n}$')
    saveFig('%i_pq_events_processed.%s' % (pltGroup, EXT))

    ax = sns.lineplot(data=lStats, x='dN', y='maxIsDelQSqrt', hue='dist', style='dist', markers=True,
                 palette=DistPalette)
    # plt.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0.8))
    plt.legend(ncol=2)
    setLegend(ax)
    plt.xscale('log')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\nicefrac{\text{events in PQ}}{\sqrt{n}}$')
    saveFig('%i_pq_events_inQueue_noInput.%s' % (pltGroup, EXT))

    ax = sns.lineplot(data=lStats, x='dN', y='maxSlSizePN', hue='dist', style='dist', markers=True,
                 palette=DistPalette)
    # plt.legend(loc='upper left', ncol=2, bbox_to_anchor=(0, 1.15))
    # plt.legend(ncol=2)
    unsetLegend(ax)
    plt.xscale('log')
    # plt.ylim((3, 10))
    plt.yscale('log', base=10)
    plt.xlabel('$n$')
    plt.ylabel(r'$\nicefrac{\text{rays in SL}}{n}$')
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

TargetNs = pd.Series(gData[gData['dist'] == 'uni']['n'].unique())


def custom_round(x):
    return TargetNs[TargetNs.sub(x).abs().idxmin()]


gData['dN'] = gData['n'].apply(lambda x: custom_round(x))

gData['tpn'] = gData['t'] / gData['n']
gData['speedup'] = 1

# calculate speedup over naive
for dist in gData['dist'].unique():
    for alg in gData['algorithm'].unique():
        for kernel in gData[gData['algorithm'] == alg]['kernel'].unique():
            maskNav = (gData['dist'] == dist) & (gData['algorithm'] == 'NaiveYao') & (gData['kernel'] == kernel)

            Ns = gData.loc[maskNav, 'n'].unique()

            maskAlg = (gData['dist'] == dist) & (gData['algorithm'] == alg) & (gData['kernel'] == kernel) & (
                gData['n'].isin(Ns))

            # gData.loc[maskAlg, 'speedup'] = gData.loc[maskNav, 't'].array / gData.loc[maskAlg, 't']

TIMELIMIT = 1800000  # ms (30 min)
Ns = gData['n'].unique()
TLpN = TIMELIMIT / Ns

######################################################################################################
# %% Runtime plots - overview plot

pltGroup = 40
ax = sns.lineplot(data=gData, x='dN', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash,
                  markers=AlgorithmMarkers,
                  palette=KernelPalette)
w, h = plt.gcf().get_size_inches()
plt.gcf().set_size_inches(1.5 * w, h)
plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.1))
setLegend(ax)
plt.xscale('log')
plt.yscale('log', base=2)
plt.xlabel(r'$n$')
plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
saveFig('%i_runtime.%s' % (pltGroup, EXT))

ax = sns.lineplot(data=gData, x='dN', y='speedup', style='kernel', hue='algorithm', dashes=KernelDash,
                  markers=KernelMarkers,
                  palette=AlgorithmPalette)
w, h = plt.gcf().get_size_inches()
plt.gcf().set_size_inches(1.5 * w, h)
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
ax = sns.lineplot(data=fData, x='dN', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash,
                  markers=AlgorithmMarkers,
                  palette=KernelPalette)
unsetLegend(ax)
plt.xscale('log')
plt.yscale('log', base=2)
plt.xlabel(r'$n$')
plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
saveFig('%i_runtime.%s' % (pltGroup, EXT))

ax = sns.lineplot(data=fData, x='dN', y='speedup', style='kernel', hue='algorithm', dashes=KernelDash,
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
    fData = gData[(gData['kernel'] == kernel) & (gData['dist'] != 'bubbles')]

    ax = sns.lineplot(data=fData, x='dN', y='tpn', style='algorithm', hue='algorithm',
                      dashes=AlgorithmDash,
                      markers=AlgorithmMarkers,
                      palette=AlgorithmPalette,
                      markersize=8)
    setLegend(ax)
    plt.xscale('log')
    plt.yscale('log', base=10)

    l, h = plt.ylim()
    plt.plot(Ns, TLpN, c='gray')
    plt.ylim(l, h)

    plt.xlabel(r'$n$')
    plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
    # plt.text(figTextX, figTextY, 'Kernel: %s' % (LegendLabels[kernel]),
    #          horizontalalignment='right',
    #          verticalalignment='bottom',
    #          transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime.%s' % (pltGroup, kernel, EXT))

    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True)

    sns.lineplot(data=fData[fData['dist'].isin(['gaussian', 'grid', 'uni'])], x='dN', y='tpn', style='algorithm', hue='dist', dashes=AlgorithmDash,
                      markers=AlgorithmMarkers,
                      palette=DistPalette,
                      markersize=8,
                      ax=axs[0])

    axs[0].set_xlabel(r'$n$')
    axs[0].set_ylabel(r'$\nicefrac{t}{n}$ [ms]')

    lH0, lL0 = (axs[0].get_legend_handles_labels())

    sns.lineplot(data=fData[fData['dist'].isin(['road', 'stars', 'circle'])], x='dN', y='tpn', style='algorithm', hue='dist', dashes=AlgorithmDash,
                      markers=AlgorithmMarkers,
                      palette=DistPalette,
                      markersize=8,
                      ax=axs[1])

    lH1, lL1 = axs[1].get_legend_handles_labels()

    unsetLegend(axs[0])
    unsetLegend(axs[1])

    lH1.insert(1, lH0[1])
    lL1.insert(1, lL0[1])

    lH1.insert(2, lH0[2])
    lL1.insert(2, lL0[2])

    lH1.insert(3, lH0[3])
    lL1.insert(3, lL0[3])

    w, h = fig.get_size_inches()
    fig.set_size_inches(2.5 * w, 1.2 * h)
    plt.legend(lH1, lL1, loc='upper left', bbox_to_anchor=(1.0, 1.05))
    setLegend(axs[1])

    plt.xscale('log')
    plt.yscale('log', base=10)

    l, h = axs[1].get_ylim()
    axs[1].plot(Ns, TLpN, c='gray')
    axs[0].plot(Ns, TLpN, c='gray')
    axs[1].set_ylim(l, h)

    axs[1].set_xlabel(r'$n$')
    # plt.ylabel(r'$\nicefrac{t}{n}$ [ms]')
    # plt.text(figTextX, figTextY, 'Kernel: %s' % (LegendLabels[kernel]),
    #          horizontalalignment='right',
    #          verticalalignment='bottom',
    #          transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime_dist.%s' % (pltGroup, kernel, EXT))

######################################################################################################
# %% Runtime plots - per kernel plots (no Naive/CGAL)
pltGroup = 55

fData = gData[(gData['algorithm'] == 'Sweepline') ^ (gData['algorithm'] == 'GridYao')]
for kernel in gData['kernel'].unique():
    # filter data by dist
    ffData = fData[fData['kernel'] == kernel]

    sns.lineplot(data=ffData, x='dN', y='tpn', style='algorithm', hue='algorithm', dashes=AlgorithmDash,
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

    sns.lineplot(data=ffData, x='dN', y='tpn', style='algorithm', hue='dist', dashes=AlgorithmDash,
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

    sns.lineplot(data=fData, x='dN', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash,
                 markers=AlgorithmMarkers,
                 palette=KernelPalette)
    plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime.%s' % (pltGroup, dist, EXT))

    for alg in fData['algorithm'].unique():
        sns.lineplot(data=fData[fData['algorithm'] == alg], x='dN', y='tpn', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
                 horizontalalignment='right',
                 verticalalignment='bottom',
                 transform=plt.gcf().transFigure)
        saveFig('%i_%s_runtime_alg_%s.%s' % (pltGroup, dist, alg, EXT))

    for kernel in fData['kernel'].unique():
        sns.lineplot(data=fData[fData['kernel'] == kernel], x='dN', y='tpn', style='algorithm', hue='kernel',
                     dashes=AlgorithmDash,
                     markers=AlgorithmMarkers, palette=KernelPalette)
        plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
                 horizontalalignment='right',
                 verticalalignment='bottom',
                 transform=plt.gcf().transFigure)
        saveFig('%i_%s_runtime_kernel_%s.%s' % (pltGroup, dist, kernel, EXT))

    sns.lineplot(data=fData[~fData['algorithm'].isin(['NaiveYao', 'CGALYao'])], x='dN', y='tpn', style='algorithm',
                 hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers, palette=KernelPalette)
    plt.text(figTextX, figTextY, 'Distribution: %s' % (dist),
             horizontalalignment='right',
             verticalalignment='bottom',
             transform=plt.gcf().transFigure)
    saveFig('%i_%s_runtime_sweepline_grid.%s' % (pltGroup, dist, EXT))

    sns.lineplot(
        data=fData[
            ~fData['algorithm'].isin(['NaiveYao', 'CGALYao']) & ~fData['kernel'].isin(['CGALExactPredExactCon'])],
        x='dN', y='tpn', style='algorithm', hue='kernel', dashes=AlgorithmDash, markers=AlgorithmMarkers,
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
plt.legend(loc='lower center', ncol=2, bbox_to_anchor=(.5, .8))
setLegend(ax)
# plt.xscale('log')
plt.yscale('log', base=10)
plt.xlabel(r'$k$')
plt.ylabel(r'$\nicefrac{t}{k}$ [ms]')
saveFig('%i_runtime_cones.%s' % (pltGroup, EXT))
