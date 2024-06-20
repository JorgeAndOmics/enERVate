from Bio import SeqIO, SeqRecord, Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

import defaults
import os
import subprocess
import concurrent.futures
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
    e_value = defaults.E_VALUE
    tmp_dir = defaults.TMP_DIR
    virus_db = defaults.VIRUS_DB
    local_db_path = defaults.VIRUS_FASTA  # Path to the local BLAST database files

def get_local_sequence(accession_id):
    """Retrieve a sequence from the local BLAST database using the accession ID."""
    db_path = Datahub.local_db_path
    try:
        with open(db_path, "r") as db_handle:
            for record in SeqIO.parse(db_handle, "fasta"):
                if record.id == accession_id:
                    return record
    except FileNotFoundError:
        logging.error(f'Local BLAST database file {db_path} not found')
        return None
    logging.error(f'Accession ID {accession_id} not found in local BLAST database')
    return None

def process_file(virus_file, key):
    file_score_list = defaultdict(list)
    try:
        with open(virus_file) as handle:
            seq_record = SeqIO.read(handle, 'fasta')

        # Construct blastn command
        logging.debug(f'Attempting to perform BLASTn:\n {seq_record}')
        blastn_cmd = [
            'blastn',
            '-query', virus_file,
            '-db', Datahub.virus_db,
            '-evalue', str(Datahub.e_value),
            '-outfmt', '5',
            '-out', os.path.join(Datahub.tmp_dir, f'tmp_result_{key}.xml')
        ]

        # Run blastn command
        result = subprocess.run(blastn_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f'Error running blastn: {result.stderr}')
            return key, None

        # Parse the BLAST results
        with open(os.path.join(Datahub.tmp_dir, f'tmp_result_{key}.xml')) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                if not blast_record.alignments:
                    logging.warning('No alignments found.')
                    continue
                for alignment in blast_record.alignments:
                    accession_id = alignment.hit_def.split()[0]  # Assuming the accession ID is the first part of the hit_def
                    for hsp in alignment.hsps:
                        logging.info(f'HSP: {hsp}')
                        file_score_list[key].append(hsp.bits)
                        logging.debug(f'Bitscore: {hsp.bits}')

                        # Retrieve the full sequence from the local BLAST database
                        full_seq_record = get_local_sequence(accession_id)
                        if not full_seq_record:
                            logging.error(f'Failed to retrieve sequence for accession ID {accession_id}')
                            continue

                        # Extract the HSP region from the full sequence
                        hsp_seq = full_seq_record.seq[hsp.sbjct_start-1:hsp.sbjct_end]

                        # Debugging log to check extracted sequence
                        logging.debug(f'Extracted HSP sequence: {hsp_seq}')

                        # Construct the output file name
                        header_parts = virus_file.split('/')[-1].replace('.fasta', '').split(';')
                        if len(header_parts) >= 5:
                            probe = header_parts[3]
                            species = header_parts[2]
                            virus_name = header_parts[1]
                            virus_family = header_parts[0]
                            n_field = header_parts[4]
                            start_location = hsp.sbjct_start
                            end_location = hsp.sbjct_end
                            strand = '+' if hsp.frame[1] > 0 else '-'
                            formatted_accession_id = full_seq_record.id.replace('.', '_')

                            output_file_name = (f"{virus_family};{virus_name};{species};{probe};{n_field};"
                                                f"st-{start_location};end-{end_location};{strand};"
                                                f"btsc-{int(hsp.bits)};{formatted_accession_id}.fasta")

                            output_dir = os.path.join('data', 'fastas', 'blastn_virus', probe, species, virus_name)
                            os.makedirs(output_dir, exist_ok=True)
                            output_file_path = os.path.join(output_dir, output_file_name)

                            # Save the HSP sequence
                            hsp_record = SeqRecord(
                                Seq(hsp_seq),
                                id=full_seq_record.id,
                                description=alignment.hit_def
                            )
                            with open(output_file_path, "w") as output_handle:
                                SeqIO.write(hsp_record, output_handle, "fasta")
                                logging.info(f'Saved HSP sequence to {output_file_path}')

        if len(file_score_list[key]) > 0:
            return key, max(file_score_list[key])
        else:
            logging.warning(f'No bitscores found for file')
            return key, None

    except (IndexError, ValueError) as e:
        logging.error(f'Error processing file {virus_file}: {e}')
        return key, None
    except FileNotFoundError as e:
        logging.error(f'File not found: {e}')
        return key, None
    except Exception as e:
        logging.exception(f'Unexpected error: {e}')
        return key, None

def virus_blaster():
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for probe in Datahub.probes:
            for species in Datahub.full_species:
                try:
                    species_dir = os.path.join(Datahub.merged_blast_dir, probe, species)
                    if not os.path.exists(species_dir):
                        logging.warning(f'No directory found for {probe} and {species}')
                        continue

                    for virus in os.listdir(species_dir):
                        virus_dir = os.path.join(species_dir, virus)
                        for file in os.listdir(virus_dir):
                            virus_file = os.path.join(virus_dir, file)
                            logging.debug(f'Processing {file}')
                            key = file.replace('.fasta', '')
                            futures.append(executor.submit(process_file, virus_file, key))

                except FileNotFoundError as e:
                    logging.warning(f'Directory not found: {e}')
                except Exception as e:
                    logging.exception(f'Unexpected error: {e}')

        file_score_dict = {}
        for future in concurrent.futures.as_completed(futures):
            key, score = future.resul
