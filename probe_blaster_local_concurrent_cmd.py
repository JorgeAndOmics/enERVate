from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

import defaults
import os
import subprocess
import concurrent.futures
import pickle
import re
from time import sleep
from urllib.error import HTTPError

import logging, coloredlogs

# Configure coloredlogs with the custom field and level styles
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.LOG_DIR, 'tblastn_log.txt')),
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
    Entrez_email: str = defaults.ENTREZ_EMAIL
    bat_species: list = defaults.BAT_SPECIES
    mammalian_species: list = defaults.MAMMAL_SPECIES
    probe_min_length: dict = defaults.PROBE_MIN_LENGTH
    expansion_switch: str = defaults.EXPANSION_SWITCH  # If 'Y', it expands the returned sequences by 'expansion_size'
    expansion_size: int = defaults.EXPANSION_SIZE
    e_value = defaults.E_VALUE
    accession_id_regex = defaults.ACCESSION_ID_REGEX
    accession_id_pattern = re.compile(accession_id_regex)
    seq_records: dict = {}  # Store the sequence records
    full_species: list = bat_species + mammalian_species
    species_db = defaults.SPECIES_DB
    full_species_taxid: dict = {}
    probes: list = defaults.PROBES
    tmp_dir = defaults.TMP_DIR
    with open(os.path.join(defaults.PICKLE_DIR, 'full_probes.pkl'), 'rb') as input_file:
        probes_dict = pickle.load(input_file)


    tblastn_results: dict = {}

def get_tax_id():
    Entrez.email = Datahub.Entrez_email
    for species_name in Datahub.full_species:
        search = Entrez.esearch(term=species_name, db='taxonomy')
        record = Entrez.read(search)
        tax_id = record['IdList'][0] if record['IdList'] else None
        Datahub.full_species_taxid[species_name] = tax_id

def perform_blast(seq_record, species, probe, virus_family, virus_name, virus_find_counter, db_path, e_value=Datahub.e_value):
    Entrez.email = Datahub.Entrez_email
    try:
        query_file = os.path.join(Datahub.tmp_dir, f'tmp_query.fasta')
        SeqIO.write(seq_record, query_file, 'fasta')

        output_file = os.path.join(Datahub.tmp_dir, f'tmp_result.xml')

        # Construct the tblastn command
        tblastn_cmd = [
            'tblastn',
            '-query', query_file,
            '-db', db_path,
            '-evalue', str(e_value),
            '-outfmt', '5',  # XML output format
            '-out', output_file
        ]

        # Run the tblastn command
        result = subprocess.run(tblastn_cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running tblastn: {result.stderr}')
            return

        # Parse the BLAST results
        with open(output_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                logging.debug(blast_record)
                if not blast_record.alignments:
                    logging.warning('No alignments found.')
                    continue
                for alignment in blast_record.alignments:
                    logging.info(f'Alignment: {alignment}')
                    for hsp in alignment.hsps:
                        logging.info(f'HSP: {hsp}')
                        hsp_length = hsp.sbjct_end - hsp.sbjct_start + 1
                        if hsp_length >= Datahub.probe_min_length[probe]:
                            accession_number = (Datahub.accession_id_pattern.search(alignment.title.split('|')[-1])).group()
                            logging.info(f'Alignment Title: {alignment.title}')
                            logging.info(f'Accession Number: {accession_number}')
                            key = (f'{virus_family};'
                                   f'{virus_name};'
                                   f'{species};'
                                   f'{probe};'
                                   f'n{virus_find_counter};'
                                   f'st-{hsp.sbjct_start};'
                                   f'end-{hsp.sbjct_end};'
                                   f's{hsp.frame[-1]};'
                                   f'{accession_number.replace(".", "_")}')

                            Datahub.tblastn_results[key] = hsp

                            # Fetch additional information using Entrez with retry logic
                            fetch_record = None
                            retries = 3
                            for attempt in range(retries):
                                try:
                                    with Entrez.efetch(db='nucleotide',
                                                       id=accession_number,
                                                       rettype='fasta',
                                                       retmode='text',
                                                       seq_start=hsp.sbjct_start,
                                                       seq_stop=hsp.sbjct_end) as handle:
                                        fetch_record = SeqIO.read(handle, 'fasta')
                                    break
                                except HTTPError as e:
                                    if attempt < retries - 1:
                                        logging.warning(f'HTTP error {e.code} encountered. Retrying... ({attempt + 1}/{retries})')
                                        sleep(2 ** attempt)  # Exponential backoff
                                    else:
                                        logging.error(f'Failed to fetch data after {retries} attempts. Error: {e}')
                                        return

                            if fetch_record:
                                Datahub.seq_records[key] = fetch_record.seq
                                expected_length = hsp.sbjct_end - hsp.sbjct_start + 1
                                actual_length = len(fetch_record.seq)
                                # Save the FASTA file
                                fasta_dir = os.path.join('data', 'fastas', 'tblastn', probe, species, virus_name)
                                if not os.path.exists(fasta_dir):
                                    os.makedirs(fasta_dir)
                                fasta_path = os.path.join(fasta_dir, f'{key}.fasta')
                                with open(fasta_path, 'w') as output_file:
                                    SeqIO.write(fetch_record, output_file, 'fasta')
                                logging.debug(f'Virus Family: {virus_family}')
                                logging.debug(f'Virus Name: {virus_name}')
                                logging.debug(f'Species: {species}')
                                logging.debug(f'Probe: {probe}')
                                logging.debug(f'Blast Number: {virus_find_counter}')
                                logging.debug(f'Alignment: {alignment}')
                                logging.debug(f'Subject Start: {hsp.sbjct_start}')
                                logging.debug(f'Subject End: {hsp.sbjct_end}')
                                logging.debug(f'Expected Length: {expected_length}, Actual Length: {actual_length}')
                                if expected_length != actual_length:
                                    logging.warning(f'Length mismatch for {key}: Expected {expected_length}, Got {actual_length}')
                                logging.debug(f'Frame: {hsp.frame}')
                                logging.info(f'FASTA file with Accession {accession_number} ({virus_name} {virus_find_counter})'
                                      f' saved at {fasta_path}.')
                                logging.debug(f'#########################################\n\n')
                                virus_find_counter += 1

        # Clean up temporary files
        os.remove(str(query_file))
        os.remove(str(output_file))
    except Exception as e:
        logging.warning(f'Warning: {str(e)}')

def fasta_blaster(db_path, e_value=0.009):
    tasks = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for species in Datahub.full_species:
            for probe in Datahub.probes:
                file_list = os.listdir(os.path.join('data', 'fastas', 'probes', probe))
                virus_find_counter = 1
                for family_virus_probe in file_list:
                    virus_family, virus_name, virus_probe = family_virus_probe.split(';')
                    for seq_record in SeqIO.parse(
                            os.path.join('data', 'fastas', 'probes', probe, family_virus_probe), 'fasta'):
                        logging.debug(f'Attempting to tBLASTn {probe} probe from {virus_name} in {species}')
                        tasks.append(
                            executor.submit(perform_blast, seq_record, species, probe, virus_family, virus_name,
                                            virus_find_counter, db_path, e_value))

        for future in concurrent.futures.as_completed(tasks):
            future.result()  # Ensure any exceptions are raised

    # Save the results to a pickle file
    with open(os.path.join('data', 'pickles', 'tblastn_results.pkl'), 'wb') as output_file:
        pickle.dump(Datahub.tblastn_results, output_file)

if __name__ == '__main__':
    # Specify the path to your local BLAST database
    local_db_path: str = Datahub.species_db
    fasta_blaster(local_db_path)
    logging.info('BLAST search completed.')
