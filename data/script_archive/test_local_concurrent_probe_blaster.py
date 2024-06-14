import os
import pickle
from Bio import Entrez, SeqIO
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import concurrent.futures


class Datahub:
    seq_records = {}  # Store the sequence records
    db_path = os.path.join('V:/databases/nucleotide_db/nucleotide_db_extracted/')  # Path to the BLAST database

def fasta_parser():
    '''
    Parse the fasta file and store the sequence records
    :return: A dictionary in Datahub containing the sequence records under their respective IDs
    '''
    # Fasta path
    full_fasta = os.path.join('fastas', 'full.fasta')

    # Read the fasta file
    for seq_record in SeqIO.parse(full_fasta, 'fasta'):
        Datahub.seq_records[seq_record.id] = seq_record.seq
        print(f'ID: {seq_record.id}')
        print(f'Description: {seq_record.description}')
        print(f'Length: {len(seq_record)}')
        print(f'Sequence: {seq_record.seq}')
        print(f'\n')

def primary_blaster_concurrent(seq=Datahub.seq_records, len_threshold=400, e_value=0.009):
    '''
    Performs a local BLAST search using the parsed contents from the full FASTA file
    :return: A FASTA file containing the sequences returned by the BLAST search
    '''

    # Perform the BLAST search
    # outfmt=5 returns the results in XML format
    # outfmt=15 returns the results in JSON format
    tblastn_cline = NcbitblastnCommandline(query=seq,
                                           db=Datahub.db_path,
                                           evalue=e_value,
                                           outfmt=5,
                                           strand='both')

    stdout, stderr = tblastn_cline()
    return stdout


def blast_parallel(num_threads, seq=Datahub.seq_records):

    # Use concurrent.futures to run BLAST in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(primary_blaster_concurrent(seq)) for id, seq in Datahub.seq_records.items()]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
        print(results)

if __name__ == '__main__':
    fasta_parser()
    blast_parallel(10)
    # print(os.listdir('V:/databases/nucleotide_db/nucleotide_db_extracted/'))


