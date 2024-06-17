import coloredlogs
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
import pandas as pd

import defaults
import os
import subprocess
import concurrent.futures
import pickle
import coloredlogs, logging
from collections import defaultdict


import logging, coloredlogs

# Configure coloredlogs with the custom field and level styles
coloredlogs.install(
    level='DEBUG',
    fmt='%(asctime)s %(hostname)s %(message)s',
    level_styles=defaults.LEVEL_STYLES,
    field_styles=defaults.FIELD_STYLES
)

class Datahub:
        bat_species = defaults.BAT_SPECIES
        mammalian_species = defaults.MAMMAL_SPECIES
        probes = defaults.PROBES
        full_species = bat_species + mammalian_species
        primary_blast_dir = os.path.join(default.FASTA_DIR, 'tblastn')
        merged_blast_dir = os.path.join(default.FASTA_DIR, 'tblastn_merged')
        virus_db = r'\\wsl$\Ubuntu\home\biouser\.ervin\virus_db_store\Viruses'


def virus_blaster():
    for probe in Datahub.probes:
        for species in Datahub.full_species:
            if not os.path.exists(os.path.join(Datahub.merged_blast_dir, probe, species)):
                logging.warning(f'No directory found for {probe} and {species}')
                continue
            else:
                for virus in os.listdir(os.path.join(Datahub.primary_blast_dir, probe, species)):

                    file_score_dict = {}

                    for file in os.listdir(os.path.join(Datahub.primary_blast_dir, probe, species, virus)):
                        virus_file = os.path.join(os.path.join(Datahub.primary_blast_dir, probe, species, virus, file))

                        key = file.split('.')[0]
                        virus_family = file.split(';')[0]
                        virus_name = file.split(';')[1]
                        species = file.split(';')[2]
                        probe = file.split(';')[3]
                        blast_number = file.split(';')[4].split('n')[1]
                        subject_start = int(file.split(';')[5].split('-')[1])
                        subject_end = int(file.split(';')[6].split('-')[1])
                        accession_id = file.split(';')[7].split('.')[0].replace("_", ".")

                        file_score_list = defaultdict(list)

                        with open(virus_file) as handle:
                            seq_record = SeqIO.read(handle, 'fasta')


                            # Construct blastn command
                            blastn_cmd = [
                                'blastn',
                                '-query', virus_file,
                                '-db', Datahub.virus_db,
                                '-evalue', str(defaults.E_VALUE),
                                '-outfmt', '5',
                                '-out', os.path.join('data', 'tmp', 'tmp_result.xml')
                            ]

                            # Run blastn command
                            result = subprocess.run(blastn_cmd, capture_output=True, text=True)

                            if result.returncode != 0:
                                logging.error(f'Error running blastn: {result.stderr}')
                                return

                            # Parse the BLAST results
                            with open(os.path.join('data', 'tmp', 'tmp_result.xml')) as result_handle:
                                blast_records = NCBIXML.parse(result_handle)
                                for blast_record in blast_records:
                                    if not blast_record.alignments:
                                        logging.warning('No alignments found.')
                                        continue
                                    for alignment in blast_record.alignments:
                                        for hsp in alignment.hsps:
                                            file_score_list[key].append(hsp.bitscore)
                        file_score_dict[key] = max(file_score_list[key])



virus_blaster()