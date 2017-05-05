from __future__ import absolute_import, print_function, division
from six.moves import range

import collections
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches


_, filename = sys.argv

records = collections.defaultdict(list)

with open(filename) as f:
    for line in f.readlines():
        form, dim, fc, time, size = line.split()
        records[(form, int(dim), fc)].append(float(time))


def select(form, dim, fc):
    return records[(form, dim, fc)]


def plot(form, dim, axes=None, left_clutter=True, right_clutter=True, top_clutter=False, bottom_clutter=False):
    if axes is None:
        axes = plt.axes()

    ticks = ['fd\_bendy', 'quadrature', 'UFLACS', 'TSFC$\dagger$', 'TSFC']

    if right_clutter:
        twinx = axes.twinx()
        twinx.get_yaxis().set_ticks([])
        l = twinx.set_ylabel('-'.join(part.capitalize() for part in form.split('_')), rotation=-90, va='bottom')
        twinx.spines['top'].set_visible(False)
        twinx.spines['left'].set_visible(False)
        twinx.spines['right'].set_visible(False)
        twinx.spines['bottom'].set_visible(False)

    if top_clutter:
        axes.set_title("{}D".format(dim))

    axes.set_xlim(xmin=0.5, xmax=len(ticks)+0.5)
    axes.get_xaxis().set_ticks(list(range(1, len(ticks)+1)))
    if bottom_clutter:
        axes.get_xaxis().set_ticklabels(ticks, rotation=40)
    else:
        axes.get_xaxis().set_ticklabels([])
    axes.yaxis.grid(True, linestyle='dotted')

    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)

    axes.set_yscale("log")
    axes.set_ylim(ymin=1e-1, ymax=1e3)
    if left_clutter:
        def myLogFormat(y, pos):
            # Find the number of decimal places required
            decimalplaces = int(np.maximum(-np.log10(y), 0))  # =0 for numbers >=1
            # Insert that number into a format string
            formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
            # Return the formatted tick label
            return formatstring.format(y)
        axes.yaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))
    else:
        axes.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ""))

    for i, fc in enumerate(['ffc-nonaffine', 'quadrature', 'uflacs', 'tsfcrepr', 'tsfc-default']):
        times = select(form, dim, fc)
        if not times:
            axes.add_patch(patches.Rectangle((i+0.5, 0.1), 1, 999.99, facecolor='gray', alpha=0.5, edgecolor="none"))
            continue

        vp = axes.violinplot([times], positions=[i+1], showmedians=True, showextrema=False)
        vp['cmedians'].set_color('black')
        for pc in vp['bodies']:
            pc.set_edgecolor('black')
            pc.set_facecolor('gray')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='10')

fig = plt.figure(figsize=(4.8, 7.4))
plot('helmholtz', 2, fig.add_subplot(421), right_clutter=False, top_clutter=True)
plot('helmholtz', 3, fig.add_subplot(422), left_clutter=False, top_clutter=True)
plot('elasticity', 2, fig.add_subplot(423), right_clutter=False)
plot('elasticity', 3, fig.add_subplot(424), left_clutter=False)
plot('hyperelasticity', 2, fig.add_subplot(425), right_clutter=False)
plot('hyperelasticity', 3, fig.add_subplot(426), left_clutter=False)
plot('holzapfel_ogden', 2, fig.add_subplot(427), right_clutter=False, bottom_clutter=True)
plot('holzapfel_ogden', 3, fig.add_subplot(428), left_clutter=False, bottom_clutter=True)
fig.tight_layout()
sp = fig.add_subplot(111, frameon=False)
sp.get_xaxis().set_visible(False)
sp.get_yaxis().set_ticks([])
sp.set_ylabel("Form compilation time [$s$]", labelpad=40)
# fig.show()
fig.savefig("example.pgf", bbox_inches='tight')
