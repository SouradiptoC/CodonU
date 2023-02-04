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
             folder_path: str = ''):
    """
    Plots ENc value against GC3 values

    :param enc_val_lst: Values of ENc
    :param gc_val_lst: Values of GC3
    :param organism_name: Name of organism (optional)
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    """
    x = list(np.linspace(0, 1, 201))
    y = [_enc(i) for i in x]
    N = len(enc_val_lst)
    color = enc_val_lst
    plt.figure(figsize=(9, 5.25))
    plt.plot(x, y, color='red', label=r"$EN_c = 2 + s + \frac{29}{s^2 + (1 - s^2)}$")
    plt.scatter(gc_val_lst, enc_val_lst, s=5, label=r"Measured $EN_c$ values", c=color, cmap='viridis')
    suptitle = r'$EN_c$ plot' if organism_name is None else f"$EN_c$ plot for {organism_name}"
    plt.suptitle(suptitle, fontsize=16)
    plt.title(f'Total genes: {N}', fontsize=12)
    plt.legend()
    plt.xlabel(r"$GC_3$ Value")
    plt.ylabel(r"$EN_c$ value")
    c_bar = plt.colorbar()
    c_bar.set_label(r'$EN_c values$')
    if save_image:
        name = 'ENc_plot.png' if organism_name is None else f"ENc_plot_{organism_name}.png"
        file_name = join(folder_path, name)
        plt.savefig(file_name, dpi=500)
    plt.show()
    plt.close()
