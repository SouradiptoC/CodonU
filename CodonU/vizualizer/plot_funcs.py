import matplotlib.pyplot as plt
import numpy as np
from os.path import join, abspath
from scipy.stats import linregress
from CodonU.file_handler.internal_comp import is_file_writeable
from CodonU.file_handler import make_dir


def _enc(x: int) -> float:
    """
    Calculates Theoretical ENC value based on Wright (1989)

    :param x: GC3 value
    :return: ENc value
    """
    return 2 + x + (29 / (x ** 2 + (1 - x) ** 2))


def _plot_enc(enc_val_lst: list, gc_val_lst: list, organism_name: None | str = None, save_image: bool = False,
              folder_path: str = '', gene_analysis: bool = True):
    """
    Plots ENc value against GC3 values

    :param enc_val_lst: Values of ENc (in range 0 to 1)
    :param gc_val_lst: Values of GC3 (in range 0 to 1)
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    """
    x = list(np.linspace(0, 1, 201))
    y = [_enc(i) for i in x]
    N = len(enc_val_lst)
    color = enc_val_lst
    plt.figure(figsize=(9, 5.25))
    plt.plot(x, y, color='red', label=r"$EN_c = 2 + s + \frac{29}{s^2 + (1 - s^2)}$", zorder=3)
    plt.scatter(gc_val_lst, enc_val_lst, s=5, label=r"Measured $EN_c$ values", c=color, cmap='viridis', alpha=0.5,
                zorder=2)
    suptitle = r'$EN_c$ plot' if organism_name is None else f"$EN_c$ plot for {organism_name}"
    plt.suptitle(suptitle, fontsize=16)
    title = f'Total genes: {N}' if gene_analysis else f'Total genome: {N}'
    plt.title(title, fontsize=12)
    plt.legend()
    plt.xlabel(r"$GC_3$ Value")
    plt.ylabel(r"$EN_c$ Value")
    plt.grid(True, linestyle=':')
    c_bar = plt.colorbar()
    c_bar.set_label(r'$EN_c values$')
    if save_image:
        make_dir(folder_path)
        name = 'ENc_plot.png' if organism_name is None else f"ENc_plot_{organism_name}.png"
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)
            print(f'Saved file can be found as {abspath(file_name)}')
    plt.show()
    plt.close()


def _plot_pr2(gc_val_lst: list, at_val_lst: list, g3_val_lst: list, a3_val_lst: list, organism_name: str | None = None,
              save_image: bool = False, folder_path: str = 'Report', gene_analysis: bool = True):
    """
    Plots A3/AT3 values against G3/GC3 values

    :param gc_val_lst: Values of GC3 (in range 0 to 1)
    :param at_val_lst: Values of AT3 (in range 0 to 1)
    :param g3_val_lst: Values of G3 (in range 0 to 1)
    :param a3_val_lst: Values of A3 (in range 0 to 1)
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    """
    N = len(gc_val_lst)
    gc_arr = np.array(gc_val_lst)
    g3_arr = np.array(g3_val_lst)
    at_arr = np.array(at_val_lst)
    a3_arr = np.array(a3_val_lst)
    x = g3_arr / gc_arr
    y = a3_arr / at_arr
    # [[x_val], [y_val]]
    ag_lst = [[], []]
    ag_count = 0
    tg_lst = [[], []]
    tg_count = 0
    tc_lst = [[], []]
    tc_count = 0
    ac_lst = [[], []]
    ac_count = 0
    no_bias = [[], []]
    no_bias_count = 0
    for x_val, y_val in zip(x, y):
        if x_val > 0.5 and y_val > 0.5:
            ag_lst[0].append(x_val)
            ag_lst[1].append(y_val)
            ag_count += 1
        elif x_val > 0.5 > y_val:
            tg_lst[0].append(x_val)
            tg_lst[1].append(y_val)
            tg_count += 1
        elif x_val < 0.5 and y_val < 0.5:
            tc_lst[0].append(x_val)
            tc_lst[1].append(y_val)
            tc_count += 1
        elif x_val < 0.5 < y_val:
            ac_lst[0].append(x_val)
            ac_lst[1].append(y_val)
            ac_count += 1
        else:
            no_bias[0].append(x_val)
            no_bias[1].append(y_val)
            no_bias_count += 1
    fig, ax = plt.subplots(figsize=(9, 5.25))
    ax.scatter(ag_lst[0], ag_lst[1], s=5, label='AG', alpha=0.5, zorder=2)
    ax.scatter(tg_lst[0], tg_lst[1], s=5, label='TG', alpha=0.5, zorder=2)
    ax.scatter(tc_lst[0], tc_lst[1], s=5, label='TC', alpha=0.5, zorder=2)
    ax.scatter(ac_lst[0], ac_lst[1], s=5, label='AC', alpha=0.5, zorder=2)
    ax.scatter(no_bias[0], no_bias[1], s=5, label=f'No Bias', alpha=0.5, zorder=3)
    ax.grid(True, linestyle=':')
    ax.legend()
    ax.text(0.6, 0.95, 'AG')
    ax.text(0.4, 0.95, 'AC')
    ax.text(0.6, 0.05, 'TG')
    ax.text(0.4, 0.05, 'TC')
    ax.text(0, 0, f'AG: {ag_count}\nTG: {tg_count}\nTC: {tc_count}\nAC: {ac_count}\nNo Bias: {no_bias_count}',
            transform=ax.transAxes)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.axvline(0.5, color='grey', zorder=1)
    plt.axhline(0.5, color='grey', zorder=1)
    plt.xlabel(r"$G_3/GC_3$ Values")
    plt.ylabel(r"$A_3/AT_3$ Values")
    suptitle = r'Parity Rule 2 plot' if organism_name is None else f"Parity Rule 2 plot for {organism_name}"
    plt.suptitle(suptitle, fontsize=16)
    title = f'Total genes: {N}' if gene_analysis else f'Total genome: {N}'
    plt.title(title, fontsize=12)
    if save_image:
        make_dir(folder_path)
        name = 'PR2_plot.png' if organism_name is None else f"PR2_plot_{organism_name}.png"
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)
            print(f'Saved file can be found as {abspath(file_name)}')
    plt.show()
    plt.close()


def _plot_neutrality(gc12_lst: list, gc3_lst: list, organism_name: None | str = None, save_image: bool = False,
                     folder_path: str = 'Report', gene_analysis: bool = True):
    """
    Plots the neutrality plot

    :param gc12_lst: The list containing values of G or C at 1 or 2 positions
    :param gc3_lst: The list containing values of G or C ar 3 positions
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    :param gene_analysis: Option if gene analysis (True) or genome analysis (False) (optional)
    """
    N = len(gc3_lst)
    if max(gc12_lst) - min(gc12_lst) > max(gc3_lst) - min(gc3_lst):
        color = gc12_lst
        label = r'$GC_{12}$ Value'
    else:
        color = gc3_lst
        label = r"$GC_3$ Value"
    slope, intercept, r, p, se = linregress(gc3_lst, gc12_lst)
    fig = plt.figure(figsize=(9, 5.25))
    ax = fig.add_subplot()
    ax.set_aspect('equal', adjustable='box')
    plt.scatter(gc3_lst, gc12_lst, s=12, c=color, cmap='viridis', alpha=0.5, label='Observed Values', zorder=1)
    x_lim = max(gc3_lst) + min(gc3_lst)
    x = np.linspace(0.0, x_lim, 201)
    y = [slope * _x + intercept for _x in x]
    _label = f"$y = {round(slope, 4)}x+{round(intercept, 4)}$" if intercept >= 0 else f"$y = {round(slope, 4)}x{round(intercept, 4)}$"
    plt.plot(x, y, color='red', label=_label, zorder=2)
    plt.grid(True, linestyle=":")
    plt.legend()
    plt.xlabel('$GC_3$ Value')
    plt.ylabel('$GC_{12}$ Value')
    c_bar = plt.colorbar()
    c_bar.set_label(label)
    suptitle = r'Neutrality Plot' if organism_name is None else f"Neutrality plot for {organism_name}"
    plt.suptitle(suptitle, fontsize=16)
    title = f'Total genes: {N}, $R^2$ value: {round(r ** 2, 4)}' if gene_analysis else f'Total genome: {N}, $R^2$ value: {round(r ** 2, 4)}'
    plt.title(title, fontsize=14)
    if save_image:
        make_dir(folder_path)
        name = 'Neutrality_plot.png' if organism_name is None else f"Neutrality_plot_{organism_name}.png"
        file_name = join(folder_path, name)
        if is_file_writeable(file_name):
            plt.savefig(file_name, dpi=500)
            print(f'Saved file can be found as {abspath(file_name)}')
    plt.show()
    plt.close()
