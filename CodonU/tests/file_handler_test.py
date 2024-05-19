from CodonU.file_handler import *
import os
import filecmp
from Bio import Entrez
from Bio import SeqIO

ACCESSION_ID = 'AP009351'
ACCESSION_ID_LST = ['CP045927.1', 'FR821777.2', 'CP000253.1', 'CP065712.1', 'AP018587.1']


def test_get_gb():
    res = get_gb(ACCESSION_ID)
    handle = Entrez.efetch(db='nucleotide', id=ACCESSION_ID, rettype='gb', retmode='text')
    record = SeqIO.read(handle, 'gb')
    assert res[ACCESSION_ID].id == record.id, "get_gb is not working"


def test_get_gbs():
    res = get_gbs(ACCESSION_ID_LST)
    handles = {acs_id: Entrez.efetch(db='nucleotide', id=acs_id, rettype='gb', retmode='text') for acs_id in
               ACCESSION_ID_LST}
    records = {acs_id: SeqIO.read(handles[acs_id], 'gb') for acs_id in ACCESSION_ID_LST}
    assert [res[acs_id].id for acs_id in ACCESSION_ID_LST] == [records[acs_id].id for acs_id in
                                                               ACCESSION_ID_LST], 'get_gbs not working'


def test_write_exome():
    # delete AP009351.1.fna
    write_exome(ACCESSION_ID)
    res = get_gb(ACCESSION_ID)
    file_name = f'{res[ACCESSION_ID].id}.fna'
    assert all([
        os.path.isfile(f'Exome/{file_name}'),
        filecmp.cmp(f'Exome/{file_name}', 'Exome/test_AP009351.1.fna')
    ]), 'write_exome not working properly'


def test_write_exomes():
    write_exomes(ACCESSION_ID_LST)
    res = get_gbs(ACCESSION_ID_LST)
    file_names = [f'{res[acs_id].id}.fna' for acs_id in ACCESSION_ID_LST]
    assert all([
        all([os.path.isfile(f'Exome/{file_name}') for file_name in file_names]),
        all([filecmp.cmp(f'Exome/{file_name}', f'Exome/test_{file_name}') for file_name in file_names])
    ])


def test_write_nucleotide():
    write_nucleotide(ACCESSION_ID)
    res = get_gb(ACCESSION_ID)
    file_name = f'{res[ACCESSION_ID].id}.ffn'
    assert all([
        os.path.isfile(f'Nuc/{file_name}'),
        filecmp.cmp(f'Nuc/{file_name}', 'Nuc/test_AP009351.1.ffn')
    ]), 'write_nucleotide not working properly'


def test_write_nucleotides():
    write_nucleotides(ACCESSION_ID_LST)
    res = get_gbs(ACCESSION_ID_LST)
    file_names = [f'{res[acs_id].id}.ffn' for acs_id in ACCESSION_ID_LST]
    assert all([
        all([os.path.isfile(f'Nuc/{file_name}') for file_name in file_names]),
        all([filecmp.cmp(f'Nuc/{file_name}', f'Nuc/test_{file_name}') for file_name in file_names])
    ])
