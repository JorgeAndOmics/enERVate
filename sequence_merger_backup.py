from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import os
import shutil
import defaults
from concurrent.futures import ProcessPoolExecutor, as_completed

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
    primary_blast_dir = os.path.join('data', 'fastas', 'tblastn')
    virus_db = r'\\wsl$\Ubuntu\home\biouser\.ervin\virus_db_store\Viruses'

def parse_fasta_headers(directory):
    """Parses the headers of FASTA files in the directory and returns a list of sequences with metadata."""
    sequences = []

    try:
        for filename in os.listdir(directory):
            if filename.endswith(".fasta"):
                filepath = os.path.join(directory, filename)
                logging.debug(f"Processing file: {filepath}")
                try:
                    for record in SeqIO.parse(filepath, "fasta"):
                        header_parts = filename.replace('.fasta', '').split(';')
                        if len(header_parts) == 8:
                            try:
                                sequences.append({
                                    "virus_family": header_parts[0],
                                    "virus_name": header_parts[1],
                                    "species": header_parts[2],
                                    "probe": header_parts[3],
                                    "n field": header_parts[4],
                                    "start": int(header_parts[5].split('-')[1]),
                                    "end": int(header_parts[6].split('-')[1]),
                                    "frame" : header_parts[7],
                                    "accession_id": header_parts[8],
                                    "sequence": str(record.seq),
                                    "fullpath": filepath
                                })
                            except (IndexError, ValueError) as e:
                                logging.error(f"Error parsing header parts for file {filepath}: {e}")
                except Exception as e:
                    logging.error(f"Error reading file {filepath}: {e}")
    except Exception as e:
        logging.error(f"Error accessing directory {directory}: {e}")

    return sequences

def merge_overlapping_sequences(sequences):
    """Merges overlapping sequences and returns a list of merged sequences with updated metadata."""
    merged_sequences_dict = {}

    try:
        # Group sequences by accession_id and frame
        sequences_by_accession_and_frame = {}
        for seq in sequences:
            key = (seq["accession_id"], seq["frame"].split('s')[-1])
            if key not in sequences_by_accession_and_frame:
                sequences_by_accession_and_frame[key] = []
            sequences_by_accession_and_frame[key].append(seq)

        for (accession_id, frame), seqs in sequences_by_accession_and_frame.items():
            seqs.sort(key=lambda x: x["start"])
            logging.debug(f"Sequences sorted by start position for accession_id: {accession_id}, frame: {frame}")

            current_seq = seqs[0]
            for next_seq in seqs[1:]:
                if next_seq["start"] <= current_seq["end"]:
                    overlap_length = current_seq["end"] - next_seq["start"] + 1
                    if overlap_length > 0:
                        current_seq["sequence"] += next_seq["sequence"][overlap_length:]
                    else:
                        current_seq["sequence"] += next_seq["sequence"]
                    current_seq["end"] = max(current_seq["end"], next_seq["end"])
                    logging.debug(f"Merging sequences: {current_seq['accession_id']} and {next_seq['accession_id']}")
                else:
                    if accession_id not in merged_sequences_dict:
                        merged_sequences_dict[accession_id] = []
                    merged_sequences_dict[accession_id].append(current_seq)
                    current_seq = next_seq

            if accession_id not in merged_sequences_dict:
                merged_sequences_dict[accession_id] = []
            merged_sequences_dict[accession_id].append(current_seq)
    except Exception as e:
        logging.error(f"Error merging sequences: {e}")

    merged_sequences = [(seq, [seq["accession_id"]]) for accession_list in merged_sequences_dict.values() for seq in accession_list]
    return merged_sequences

def verify_merged_sequences(merged_sequences):
    """Verifies the merged sequences for correct length."""
    for seq, accession_ids in merged_sequences:
        expected_length = seq['end'] - seq['start'] + 1
        if len(seq['sequence']) != expected_length:
            logging.warning(f"Sequence length mismatch for {accession_ids[0]}: "
                            f"expected {expected_length}, got {len(seq['sequence'])}")
            return False
    return True

def write_merged_sequences(original_dir, target_dir, merged_sequences, all_sequences):
    """Writes the merged sequences to new FASTA files with a consistent naming convention."""
    try:
        os.makedirs(target_dir, exist_ok=True)
        logging.debug(f"Created target directory: {target_dir}")
    except Exception as e:
        logging.error(f"Error creating directory {target_dir}: {e}")
        return

    merged_accession_ids = set()

    # Write merged sequences
    for seq, accession_ids in merged_sequences:
        try:
            expected_length = seq['end'] - seq['start'] + 1
            if len(seq['sequence']) != expected_length:
                logging.warning(f"Sequence length mismatch for {accession_ids[0]}: "
                                f"expected {expected_length}, got {len(seq['sequence'])}")
                continue

            header = f"{seq['virus_family']};{seq['virus_name']};{seq['species']};{seq['probe']};n();st-{seq['start']};end-{seq['end']};{accession_ids[0]}.fasta"
            record = SeqRecord(Seq(seq["sequence"]), id=header, description="")
            output_filepath = os.path.join(target_dir, header)
            with open(output_filepath, "w") as output_file:
                SeqIO.write(record, output_file, "fasta")
            merged_accession_ids.update(accession_ids)
            logging.info(f"Wrote merged sequence to {output_filepath}")
        except Exception as e:
            logging.error(f"Error writing merged sequence: {e}")

    # Copy unmerged sequences
    for seq in all_sequences:
        if seq['accession_id'] not in merged_accession_ids:
            try:
                shutil.copy(seq['fullpath'], target_dir)
                logging.debug(f"\n Copied unmerged file \n"
                             f"Family: {seq['virus_family']} \n"
                             f"Virus: {seq['virus_name']} \n"
                             f"Host: {seq['species']} \n"
                             f"Probe: {seq['probe']} \n"
                             f"Start: {seq['start']} \n"
                             f"End: {seq['end']} \n"
                             f"Frame: {seq['frame']} \n"
                             f"Accession ID: {seq['accession_id']} \n")
            except Exception as e:
                logging.error(f"Error copying file {seq['fullpath']} to {target_dir}: {e}")

def process_directory(args):
    original_dir, target_base_dir = args
    try:
        # Creating the new target directory structure
        relative_path = os.path.relpath(original_dir, start=os.path.dirname(original_dir))
        target_dir = os.path.join(target_base_dir, relative_path)

        logging.debug(f"Starting processing for directory {original_dir}")
        sequences = parse_fasta_headers(original_dir)
        merged_sequences = merge_overlapping_sequences(sequences)

        # Verify merged sequences
        if not verify_merged_sequences(merged_sequences):
            logging.error(f"Verification failed for directory {original_dir}")
            return

        write_merged_sequences(original_dir, target_dir, merged_sequences, sequences)
        logging.info(f"Finished processing for directory {original_dir}")
    except Exception as e:
        logging.error(f"Error in main processing: {e}")

if __name__ == '__main__':
    directories_to_process = []
    for probes in Datahub.probes:
        for species in Datahub.full_species:
            if not os.path.exists(os.path.join(Datahub.primary_blast_dir, probes, species)):
                logging.warning(f'No directory found for {probes} and {species} -> Skipping...')
                continue
            for virus in os.listdir(os.path.join(Datahub.primary_blast_dir, probes, species)):
                original_directory = os.path.join(Datahub.primary_blast_dir, probes, species, virus)
                target_base_directory = os.path.join('data', 'fastas', 'tblastn_merged', probes, species)
                if not os.path.exists(target_base_directory):
                    logging.warning(f'No target directory found for {probes} and {species} -> Generating...')
                    os.makedirs(target_base_directory)
                directories_to_process.append((original_directory, target_base_directory))

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_directory, args) for args in directories_to_process]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                logging.error(f'Error in processing: {exc}')
