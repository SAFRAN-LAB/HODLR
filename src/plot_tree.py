import numpy as np
import matplotlib.patches as patches
import pylab as pl

# Setting plot parameters for nicer plots:
pl.rcParams['figure.figsize']  = 15, 10
pl.rcParams['figure.dpi']      = 500
pl.rcParams['image.cmap']      = 'jet'
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 30
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'

pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad']  = 8
pl.rcParams['xtick.minor.pad']  = 8
pl.rcParams['xtick.color']      = 'k'
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'

pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad']  = 8
pl.rcParams['ytick.minor.pad']  = 8
pl.rcParams['ytick.color']      = 'k'
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'

class Node(object):
    def __init__(self, cx, cy, rx, ry):
        self.cx = cx
        self.cy = cy
        self.rx = rx
        self.ry = ry

def return_HODLR_tree(N_levels):
    # Adding root:
    tree = [[Node(0.5, 0.5, 0.5, 0.5)]]

    for i in range(1, N_levels + 1):
        level = []
        for j in range(2**(i-1)):
            parent_box = tree[-1][j]
            child_rx = parent_box.rx / 2
            child_ry = parent_box.ry / 2

            child_cx1 = parent_box.cx - child_rx
            child_cy1 = parent_box.cy + child_ry

            child_cx2 = parent_box.cx + child_rx
            child_cy2 = parent_box.cy - child_ry

            level.append(Node(child_cx1, child_cy1, child_rx, child_ry))
            level.append(Node(child_cx2, child_cy2, child_rx, child_ry))

        tree.append(level)

    return tree

def extract_centers_radii(tree):
    cx = []
    cy = []

    rx = []
    ry = []

    for i, level in enumerate(tree):
        for j in range(2**i):
            cx.append(level[j].cx)
            cy.append(level[j].cy)
            rx.append(level[j].rx)
            ry.append(level[j].ry)

    cx = np.array(cx)
    cy = np.array(cy)
    rx = np.array(rx)
    ry = np.array(ry)

    return cx, cy, rx, ry

def plot_graph(cx, cy, rx, ry, rank):

    fig = pl.figure()

    x_min = (cx - rx).min()
    y_min = (cy - ry).min()
    x_max = (cx + rx).max()
    y_max = (cy + ry).max()

    # For some reason the plot functionality wants
    # things in the range of [0, 1]
    # Mapping values to [0, 1]:
    cx = (cx - x_min) / (x_max - x_min)
    rx = rx / (x_max - x_min)

    cy = (cy - y_min) / (y_max - y_min)
    ry = ry / (y_max - y_min)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    ax.set_aspect('equal')

    # Asserting that they are all of the same length:
    if(not (len(cx) == len(cy) == len(rx) == len(ry))):
        raise AssertionError('Elements in array are not the same!!')

    for i, cx_node in enumerate(cx):
        if(rx[i] == np.min(rx) and cx[i] == 1 - cy[i]):
            ax.add_patch(patches.Rectangle((cx_node - rx[i], cy[i] - ry[i]),
                                           2 * rx[i], 2 * ry[i], linewidth = 0.1, facecolor = 'red',
                                           edgecolor = 'black'))

        else:
            ax.text(cx_node + rx[i] / 2, cy[i] + ry[i] / 2, '%02d'%(rank[i]),
                    fontsize = 30 / (1 + int(np.sqrt(i))))
            ax.text(cx_node - rx[i] / 2, cy[i] - ry[i] / 2, '%02d'%(rank[i]),
                    fontsize = 30 / (1 + int(np.sqrt(i))))

            if(rank[i] > 0):
                if(np.min(rank) == np.max(rank)):
                    intensity = 1
                else:
                    intensity = 0.2 + 0.8 * ((rank[i] - np.min(rank)) / (np.max(rank) - np.min(rank)))

                ax.add_patch(patches.Rectangle((cx_node, cy[i]),
                                               rx[i], ry[i], facecolor = 'green',
                                               edgecolor = 'black', linewidth = 0.1,
                                               alpha = intensity))
                ax.add_patch(patches.Rectangle((cx_node, cy[i]),
                                               -rx[i], -ry[i], facecolor = 'green',
                                               edgecolor = 'black', linewidth = 0.1,
                                               alpha = intensity))
            else:
                ax.add_patch(patches.Rectangle((cx_node, cy[i]),
                                               rx[i], ry[i], facecolor = 'white',
                                               edgecolor = 'black', linewidth = 0.1))
                ax.add_patch(patches.Rectangle((cx_node, cy[i]),
                                               -rx[i], -ry[i], facecolor = 'white',
                                               edgecolor = 'black', linewidth = 0.1))

import sys
imgname = str(sys.argv[1])
rank    = np.loadtxt("rank.txt")
levels  = int(rank[0])
cx, cy, rx, ry = extract_centers_radii(return_HODLR_tree(levels))
plot_graph(cx, cy, rx, ry, rank[1:])
pl.savefig(imgname, bbox_inches = 'tight')
