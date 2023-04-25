from Bio.Phylo import draw, read
from os.path import join, abspath
import matplotlib.pyplot as plt
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_writeable


def plot_phy_dnd(handle: str, title: str = 'Phylogenetic Tree', save_image: bool = False, folder_path: str = 'Report'):
    """
    Plots phylogenetic tree from dnd file
    :param handle: Handle to the dnd file
    :param title: Title of the plot
    :param save_image: Options for saving the image (optional)
    :param folder_path: Folder path where image should be saved (optional)
    """
    tree = read(handle, 'newick')
    tree = tree.as_phyloxml()
    fig, ax = plt.subplots()
    fig.suptitle(title)
    draw(tree, axes=ax)
    if save_image:
        make_dir(folder_path)
        file_name = '_'.join(title.split()) + '_dnd.png'
        file_path = join(folder_path, file_name)
        if is_file_writeable(file_path):
            fig.savefig(file_path, dpi=500)
            print(f'Saved file can be found as {abspath(file_path)}')
    plt.close(fig)
