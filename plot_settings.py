import matplotlib
import numpy as np

def settings(twocolumn=False):

    if twocolumn:
        fig_width_pt = 255.76535  # Get this from LaTeX using \showthe\columnwidth
        fontzise = 10
    else:
        fig_width_pt = 426.79134  # Get this from LaTeX using \showthe\columnwidth
        fontzise = 12

    inches_per_pt = 1.0 / 72.27  # Convert pt to inch
    golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean  # height in inches
    fig_size = [fig_width, fig_height]

    params = {'backend': 'pdf',
          'axes.labelsize': 8,
          'legend.fontsize': 8,
          'font.size': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'axes.linewidth': 0.5,
          'lines.linewidth': 0.5,
          'legend.frameon': False,
          'text.usetex': True,
          'ps.usedistiller': False,
          'figure.figsize': fig_size,
          'font.family': 'Computer Modern',
          'font.fantasy': ['Comic Sans MS'],
          }
    matplotlib.rcParams.update(params)
    return ()