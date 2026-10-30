import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
    'grid.linewidth': 0.6, 'axes.axisbelow': True, 'legend.frameon': False, 'lines.linewidth': 2,
    'savefig.dpi': 300, 'savefig.bbox': 'tight', 'figure.dpi': 150})


def panel(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom', ha='left', color=INK)
