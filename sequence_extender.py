from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

import defaults
import os
import subprocess
import pickle
import logging

# Configure coloredlogs with the custom field and level styles
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.LOG_DIR, 'extending_log.txt')),
        logging.StreamHandler()
             ]
                    )

coloredlogs.install(
    level='DEBUG',
    fmt='%(asctime)s - %(message)s',
    level_styles=defaults.LEVEL_STYLES,
    field_styles=defaults.FIELD_STYLES
)

class Datahub:
    expansion_size = EXPANSION_SIZE


def seq_expander(expansion_size=Datahub.expansion_size):
    Entrez.email = 'jgonzlez@tcd.ie'
    for dirpaths, dirnames, filenames in os.walk(os.path.join('data', 'fastas', 'tblastn')):
        for dir_count, file in enumerate(filenames):
            try:
                logging.info(f'Expanding {file} by {expansion_size} nucleotides')
                # Fetch the nucleotide record
                virus_family = file.split(';')[0]
                virus_name = file.split(';')[1]
                species = file.split(';')[2]
                probe = file.split(';')[3]
                blast_number = file.split(';')[4].split('n')[1]
                subject_start = int(file.split(';')[5].split('-')[1])
                subject_end = int(file.split(';')[6].split('-')[1])
                reading_frame = file.split(';')[7].split('.')[1]
                accession_id = file.split(';')[8].split('.')[0].replace("_", ".")

                logging.debug(f'Virus Family: {virus_family}')
                logging.debug(f'Virus Name: {virus_name}')
                logging.debug(f'Species: {species}')
                logging.debug(f'Probe: {probe}')
                logging.debug(f'Blast Number: {blast_number}/{dir_count}')
                logging.debug(f'Subject Start: {subject_start}')
                logging.debug(f'Subject End: {subject_end}')
                logging.debug(f'Reading Frame: {reading_frame}')
                logging.debug(f'Accession ID: {accession_id}')

                with Entrez.efetch(db="nucleotide",
                                   id=accession_id,
                                   rettype="gb",
                                   retmode="text",
                                   seq_start=max(1, subject_start - expansion_size),
                                   seq_stop=subject_end + expansion_size) as handle:
                    record = SeqIO.read(handle, "genbank")

                logging.debug(f'Original Sequence Length: {subject_end - subject_start + 1}')
                logging.debug(f'Expanded Sequence Length: {len(record.seq)}')

                fasta_path = os.path.join('data',
                                          'fastas',
                                          'tblastn_expanded',
                                          probe,
                                          species,
                                          virus_name)

                os.makedirs(fasta_path, exist_ok=True)

                fasta_file = os.path.join(fasta_path, f'{file}')

                with open(fasta_file, 'w') as output_file:
                    SeqIO.write(record, output_file, 'fasta')

            except Exception as e:
                logging.error(f"An error occurred while processing {file}: {e}", exc_info=True)


if __name__ == '__main__':
    try:
        seq_expander()
    except Exception as e:
        logging.critical(f"A critical error occurred in the main execution: {e}", exc_info=True)
