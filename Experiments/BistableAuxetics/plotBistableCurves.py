import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

LINE_COLOR = '#1F5EAA'
INK_COLOR = '#404040'


def load_result(prefix):
    return loadmat(prefix + '.mat', squeeze_me=True, struct_as_record=False)['r']


def plot_panel(ax, strain, values, ylabel, title):
    ax.plot(strain, values, '-o', color=LINE_COLOR, linewidth=2, markersize=6)
    ax.set_xlim(strain.min(), strain.max())
    ax.axhline(0, color=INK_COLOR, linewidth=0.8)
    ax.set_xlabel('Deformación impuesta ε', color=INK_COLOR)
    ax.set_ylabel(ylabel, color=INK_COLOR)
    ax.set_title(title, color=INK_COLOR, fontsize=12)
    ax.grid(alpha=0.2)
    ax.tick_params(colors=INK_COLOR)
    ax.spines[['top', 'right']].set_visible(False)


def plot_energy_and_force(result, file_name):
    fig, (ax_energy, ax_force) = plt.subplots(1, 2, figsize=(11, 4.4))
    plot_panel(ax_energy, result.strain, result.energy,
               'Energía E / μ [mm²]', 'Energía de deformación')
    plot_panel(ax_force, result.strain, result.force,
               'dE/dε / μ [mm²]', 'Fuerza conjugada de ε')
    fig.tight_layout()
    fig.savefig(file_name, dpi=150)


def main(prefix):
    plot_energy_and_force(load_result(prefix), prefix + '_curves.png')


if __name__ == '__main__':
    main(sys.argv[1])
