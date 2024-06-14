import os
import pickle
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import concurrent.futures
from tqdm import tqdm


class Datahub:
    seq_records = {}  # Store the sequence records

def fasta_parser():
    '''
    Parse the fasta file and store the sequence records
    :return: A dictionary in Datahub containing the sequence records under their respective IDs
    '''
    # Fasta path
    full_fasta = os.path.join('..', 'fastas', 'full.fasta')

    # Read the fasta file
    for seq_record in SeqIO.parse(full_fasta, 'fasta'):
        Datahub.seq_records[seq_record.id] = seq_record.seq
        print(f'ID: {seq_record.id}')
        print(f'Description: {seq_record.description}')
        print(f'Length: {len(seq_record)}')
        print(f'Sequence: {seq_record.seq}')
        print(f'\n')

def primary_blaster_concurrent(seq, len_threshold=400, e_value=0.009):
    '''
    Perform a BLAST search using the parsed contents from the full FASTA file
    :return: A FASTA file containing the sequences returned by the BLAST search
    '''
    Entrez.email = 'jorgegonzalezvet@gmail.com'
    blast_sequences = {}

    # Perform the BLAST search
    result_handle = NCBIWWW.qblast('tblastn',
                                   'nt',
                                   seq,
                                   expect=e_value
                                   )

    # Parse the BLAST results
    blast_results = {}
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if alignment.title not in blast_results and alignment.length > len_threshold:
                    blast_results[alignment.title] = {
                        'title': alignment.title,
                        'length': alignment.length,
                        'score': hsp.score,
                        'e_value': hsp.expect,
                        'identities': hsp.identities,
                        'gaps': hsp.gaps,
                        'query': hsp.query,
                        'match': hsp.match,
                        'sbjct': hsp.sbjct
                    }

    # Retrieve the sequences returned by the BLAST search in FASTA format
    for title, result in blast_results.items():
        # Get the accession number from the title
        accession_number = title.split('|')[3]
        # Fetch the sequence from NCBI
        print(f'Fetching sequence for {accession_number}')
        handle = Entrez.efetch(db='nucleotide', id=accession_number, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, 'fasta')
        handle.close()
        blast_sequences['title'] = record.seq

    # Concatenate and save the sequences in blast_sequences to a FASTA file
    with open(os.path.join('..', 'fastas', 'test_concurrent_primary_blast.fasta'), 'w') as output_file:
        print(f'Saving the sequences to {output_file}')
        for title, sequence in blast_sequences.items():
            seq_record = SeqRecord(sequence, id=title, description='')
            SeqIO.write(seq_record, output_file, 'fasta')



if __name__ == '__main__':
    fasta_parser()

    # Using ThreadPoolExecutor to parallelize
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        print(f'Executing BLAST search...')
        results = list(tqdm(executor.map(primary_blaster_concurrent, Datahub.seq_records.values()),
                            total=len(Datahub.seq_records.values())))



