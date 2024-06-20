from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

import defaults
import pandas as pd
import os
import subprocess
import concurrent.futures
import pickle
from collections import defaultdict


import logging, coloredlogs

# Configure coloredlogs with the custom field and level styles
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.LOG_DIR, 'virus_blasting_log.txt'), mode='w'),
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
    bat_species = defaults.BAT_SPECIES
    mammalian_species = defaults.MAMMAL_SPECIES
    full_species = bat_species + mammalian_species
    probes = defaults.PROBES
    merged_blast_dir = os.path.join(defaults.FASTA_DIR, 'tblastn_merged')
    primary_blast_dir = os.path.join(defaults.FASTA_DIR, 'tblastn')
    tmp_dir = defaults.TMP_DIR
    virus_db = defaults.VIRUS_DB

def virus_blaster():
    for probe in Datahub.probes:
        for species in Datahub.full_species:
            try:
                species_dir = os.path.join(Datahub.merged_blast_dir, probe, species)
                if not os.path.exists(species_dir):
                    logging.warning(f'No directory found for {probe} and {species}')
                    continue

                for virus in os.listdir(species_dir):
                    virus_dir = os.path.join(species_dir, virus)
                    file_score_dict = {}

                    for file in os.listdir(virus_dir):
                        virus_file = os.path.join(virus_dir, file)
                        logging.debug(f'Processing {file}')
                        try:
                            key = file.replace('.fasta', '')
                            parts = file.split(';')
                            virus_family = parts[0]
                            virus_name = parts[1]
                            species = parts[2]
                            probe = parts[3]
                            blast_number = parts[4].split('n')[1]
                            subject_start = int(parts[5].split('-')[1])
                            subject_end = int(parts[6].split('-')[1])
                            strand = parts[7]
                            accession_id = parts[8].split('.')[0].replace("_", ".")

                            file_score_list = defaultdict(list)

                            with open(virus_file) as handle:
                                seq_record = SeqIO.read(handle, 'fasta')

                            # Construct blastn command
                            logging.debug(f'Attempting to perform BLASTn:\n {seq_record}')
                            blastn_cmd = [
                                'blastn',
                                '-query', virus_file,
                                '-db', Datahub.virus_db,
                                '-evalue', str(0.1),
                                '-outfmt', '5',
                                '-out', os.path.join(Datahub.tmp_dir, 'tmp_result.xml')
                            ]

                            # Run blastn command
                            result = subprocess.run(blastn_cmd, capture_output=True, text=True)
                            if result.returncode != 0:
                                logging.error(f'Error running blastn: {result.stderr}')
                                continue

                            # Parse the BLAST results
                            with open(os.path.join(Datahub.tmp_dir, 'tmp_result.xml')) as result_handle:
                                blast_records = NCBIXML.parse(result_handle)
                                for blast_record in blast_records:
                                    if not blast_record.alignments:
                                        logging.warning('No alignments found.')
                                        continue
                                    for alignment in blast_record.alignments:
                                        for hsp in alignment.hsps:
                                            file_score_list[key].append(hsp.bits)
                                            logging.debug(f'Bitscore: {hsp.bits}')

                            if len(file_score_list[key]) > 0:
                                file_score_dict[key] = max(file_score_list[key])
                                logging.info(f'Bitscores found for file: {file_score_dict[key]}')
                            else:
                                logging.warning(f'No bitscores found for file')

                        except (IndexError, ValueError) as e:
                            logging.error(f'Error processing file {file}: {e}')
                        except FileNotFoundError as e:
                            logging.error(f'File not found: {e}')
                        except Exception as e:
                            logging.exception(f'Unexpected error: {e}')

            except FileNotFoundError as e:
                logging.error(f'Directory not found: {e}')
            except Exception as e:
                logging.exception(f'Unexpected error: {e}')

if __name__ == "__main__":
    try:
        virus_blaster()
    except Exception as e:
        logging.critical(f'Critical error in virus_blaster: {e}')