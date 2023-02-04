import matplotlib.pyplot as plt
import numpy as np
from os.path import join


def _enc(x: int) -> float:
    """
    Calculates Theoretical ENC value based on Wright (1989)

    :param x: GC3 value
    :return: ENc value
    """
    return 2 + x + (29 / (x ** 2 + (1 - x) ** 2))


def plot_enc(enc_val_lst: list, gc_val_lst: list, organism_name: None | str = None, save_image: bool = False,
             folder_path: str = '', gene_analysis: bool = True):
    """
    Plots ENc value against GC3 values

    :param enc_val_lst: Values of ENc
    :param gc_val_lst: Values of GC3
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
    plt.plot(x, y, color='red', label=r"$EN_c = 2 + s + \frac{29}{s^2 + (1 - s^2)}$")
    plt.scatter(gc_val_lst, enc_val_lst, s=5, label=r"Measured $EN_c$ values", c=color, cmap='viridis', alpha=0.5)
    suptitle = r'$EN_c$ plot' if organism_name is None else f"$EN_c$ plot for {organism_name}"
    plt.suptitle(suptitle, fontsize=16)
    title = f'Total genes: {N}' if gene_analysis else f'Total genome: {N}'
    plt.title(title, fontsize=12)
    plt.legend()
    plt.xlabel(r"$GC_3$ Value")
    plt.ylabel(r"$EN_c$ Value")
    c_bar = plt.colorbar()
    c_bar.set_label(r'$EN_c values$')
    if save_image:
        name = 'ENc_plot.png' if organism_name is None else f"ENc_plot_{organism_name}.png"
        file_name = join(folder_path, name)
        plt.savefig(file_name, dpi=500)
    plt.show()
    plt.close()


def plot_pr2(gc_val_lst: list, at_val_lst: list, g3_val_lst: list, a3_val_lst: list, organism_name: str | None = None,
             save_image: bool = False, folder_path: str = '', gene_analysis: bool = True):
    """
    Plots A3/AT3 values against G3/GC3 values

    :param gc_val_lst: Values of GC3
    :param at_val_lst: Values of AT3
    :param g3_val_lst: Values of G3
    :param a3_val_lst: Values of A3
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
    plt.figure(figsize=(9, 5.25))
    plt.scatter(ag_lst[0], ag_lst[1], s=5, label='AG', alpha=0.5, zorder=2)
    plt.scatter(tg_lst[0], tg_lst[1], s=5, label='TG', alpha=0.5, zorder=2)
    plt.scatter(tc_lst[0], tc_lst[1], s=5, label='TC', alpha=0.5, zorder=2)
    plt.scatter(ac_lst[0], ac_lst[1], s=5, label='AC', alpha=0.5, zorder=2)
    plt.scatter(no_bias[0], no_bias[1], s=5, label=f'No Bias', alpha=0.5, zorder=3)
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.text(0.6, 0.95, 'AG')
    plt.text(0.4, 0.95, 'AC')
    plt.text(0.6, 0.05, 'TG')
    plt.text(0.4, 0.05, 'TC')
    plt.text(0.87, 0.04, f'AG: {ag_count}\nTG: {tg_count}\nTC: {tc_count}\nAC: {ac_count}\nNo Bias: {no_bias_count}',
             bbox={'facecolor': 'white', 'alpha': 0.5, 'boxstyle': 'round', 'pad': 0.5})
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.axvline(0.5, color='grey', zorder=1)
    plt.axhline(0.5, color='grey', zorder=1)
    plt.xlabel(r"$G_3/GC_3$ Values")
    plt.ylabel(r"$A_3/AT_3$ Values")
    suptitle = r'Pairty Rule 2 plot' if organism_name is None else f"Pairty Rule 2 plot for {organism_name}"
    plt.suptitle(suptitle, fontsize=16)
    title = f'Total genes: {N}' if gene_analysis else f'Total genome: {N}'
    plt.title(title, fontsize=12)
    if save_image:
        name = 'PR2_plot.png' if organism_name is None else f"PR2_plot_{organism_name}.png"
        file_name = join(folder_path, name)
        plt.savefig(file_name, dpi=500)
    plt.show()
    plt.close()
