from Bio.Align.Applications import ClustalwCommandline
from CodonU.file_handler import make_dir
from CodonU.file_handler.internal_comp import is_file_empty, is_file
from os.path import join, abspath
import subprocess


def phy_clustal_w(bin_path: str, handle: str, res_folder_path: str = 'Report'):
    """
    Makes the multiple sequence alignment with ClustalW. For details visit http://www.clustal.org/clustal2

    :param bin_path: Path to the binary file of ClustalW
    :param handle: Handle to the file, or the filename as a string
    :param res_folder_path: The path of folder where the file will be saved
    """
    make_dir(res_folder_path)
    identifier = handle.split('/')[-1].split('.')[0]
    report_file_name = f'{identifier}_aligned_w.nex'
    report_file_path = join(res_folder_path, report_file_name)
    if not is_file(report_file_path) or is_file_empty(report_file_path):
        clustalW_cline = ClustalwCommandline(bin_path, infile=handle, outfile=report_file_path, output='NEXUS',
                                             outputtree='nexus')
        cmd = clustalW_cline.__str__()
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(f'The alignment file can be can be found at: {abspath(report_file_path)}')
