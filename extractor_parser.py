import defaults
import pandas as pd
import numpy as np
import os
import pickle
from collections import defaultdict
from Bio import Entrez, SeqIO



class Datahub:
    bat_species = defaults.BAT_SPECIES
    mammalian_species = defaults.MAMMAL_SPECIES
    full_genomes = bat_species + mammalian_species
    full_probes = {}  # Store the probe genes and their respective IDs

def table_parser():
    '''
    Parse the CSV files containing the probe genes and their respective IDs
    :input: Takes in the CSV files containing the probe genes and their respective IDs in the tables directory
    :return: A dictionary containing the probe genes {family;name;probe} and their respective IDs
    '''
    full_probes = {}
    gag_probes = {}
    pol_probes = {}
    env_probes = {}
    vif_probes = {}
    n_probes = {}
    p_probes = {}
    g_probes = {}
    l_probes = {}
    x_probes = {}
    m_probes = {}

    # Read the CSV file
    table_gpe = pd.read_csv(os.path.join('data', 'tables', 'table_gag_pol_env.csv'))
    table_vif = pd.read_csv(os.path.join('data', 'tables', 'table_vif.csv'))
    table_borna = pd.read_csv(os.path.join('data', 'tables', 'table_borna.csv'))

    for index, row in table_gpe.iterrows():
        family = row['Family']
        abbreviation = row['Abbreviation']
        name = row['Name']
        pol = row['POL']
        print(f'Parsed {family};{name} POL ID: {pol}')
        env = row['ENV']
        print(f'Parsed {family};{name} ENV ID: {env}')
        gag = row['GAG']
        print(f'Parsed {family};{name} GAG ID: {gag}')

        pol_probes[f'{family};{name};pol'] = pol
        env_probes[f'{family};{name};env'] = env
        gag_probes[f'{family};{name};gag'] = gag

    for index, row in table_vif.iterrows():
        family = row['Family']
        abbreviation = row['Abbreviation']
        name = row['Name']
        vif = row['VIF']
        print(f'Parsed {family};{name} VIF ID: {vif}')

        vif_probes[f'{family};{name};vif'] = vif

    for index, row in table_borna.iterrows():
        family = row['Family']
        abbreviation = row['Abbreviation']
        name = row['Name']
        n_protein = row['N']
        print(f'Parsed {family};{name} N Protein ID: {n_protein}')
        p_protein = row['P']
        print(f'Parsed {family};{name} P Protein ID: {p_protein}')
        g_protein = row['G']
        print(f'Parsed {family};{name} G Protein ID: {g_protein}')
        l_protein = row['L']
        print(f'Parsed {family};{name} L Protein ID: {l_protein}')
        x_protein = row['X']
        print(f'Parsed {family};{name} X Protein ID: {x_protein}')
        m_protein = row['M']
        print(f'Parsed {family};{name} M Protein ID: {m_protein}')


        n_probes[f'{family};{name};n_protein'] = n_protein
        p_probes[f'{family};{name};p_protein'] = p_protein
        g_probes[f'{family};{name};g_protein'] = g_protein
        l_probes[f'{family};{name};l_protein'] = l_protein
        x_probes[f'{family};{name};x_protein'] = x_protein
        m_probes[f'{family};{name};m_protein'] = m_protein


    full_probes['pol'] = pol_probes
    full_probes['env'] = env_probes
    full_probes['gag'] = gag_probes
    full_probes['vif'] = vif_probes
    full_probes['n_protein'] = n_probes
    full_probes['p_protein'] = p_probes
    full_probes['g_protein'] = g_probes
    full_probes['l_protein'] = l_probes
    full_probes['x_protein'] = x_probes
    full_probes['m_protein'] = m_probes

    def dict_pickler():
        with open(os.path.join('data', 'pickles', 'full_probes.pkl'), 'wb') as output_file:
            pickle.dump(full_probes, output_file)
        print('Full probes dictionary pickled.')

    dict_pickler()

    return full_probes


def download_genes(dicts) -> 'Set of probe genes':
    '''
    Download the sequences for the probe genes
    :input: A dictionary containing the probe genes
    :return: A set of probe genes
    '''
    entrez_email = 'jorgegonzalezvet@gmail.com'
    for gene_type, gene_info in dicts.items():
        for gene_key, accession_number in gene_info.items():
            print(f'Downloading sequences for {gene_key} -> {accession_number}')

            # Creating the directory if it doesn't exist
            directory_path = os.path.join('data', 'fastas', 'probes', gene_type)
            if not os.path.exists(directory_path):
                os.makedirs(directory_path)

            try:
                # Fetching the record from NCBI Protein database
                handle = Entrez.efetch(db="protein", id=accession_number, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                handle.close()

                if record:
                    # Extracting relevant information from the record
                    description = record.description
                    sequence_length = len(record.seq)
                    print(f'Record: {record}')


                    # Save the sequences to a file
                    file_path = os.path.join(directory_path, f'{gene_key}.fasta')
                    with open(file_path, "w") as output_file:
                        output_file.write(record.format("fasta"))
                    print(f'Sequences for {gene_key} downloaded and saved.')

                else:
                    print(f"No record found for {gene_key}.")
                    with open(os.path.join('data', 'logs', 'missing_records.txt'), "a") as output_file:
                        output_file.write(f'> No record found for {" ".join(gene_key.split("-"))}.\n')

            except Exception as e:
                print(f'An error occurred while downloading {gene_key}: {str(e)}')
                with open(os.path.join('data', 'logs', 'missing_records.txt'), "a") as output_file:
                    output_file.write(f'> No record found for {" ".join(gene_key.split("-"))}.\n')

    print("Gene sequences have been downloaded and saved.")


def concatenate_fastas() -> 'FASTA':
    """
    INCOMPLETE Concatenate multiple FASTA files from a directory into a single FASTA file.

    Parameters:
    input_dir (str): Directory containing the input FASTA files.
    output_file (str): Path to the output FASTA file.
    """
    # Get a list of FASTA files in the specified directory
    fasta_files: list = []
    for root, dirs, files in os.walk(os.path.join('data', 'fastas', 'probes')):
        for file in files:
            if file.endswith('.fasta') or file.endswith('.fa'):
                fasta_files.append(file)


    # Open the output file in write mode
    probes: list = ['gag', 'pol', 'env']
    output_file = os.path.join('data', 'fastas')
    with open(os.path.join('data', 'fastas', 'full.fasta'), 'w') as outfile:
        # Iterate over the list of FASTA files
        for probe in probes:
            # Open the current FASTA file in read mode
            for file in os.listdir(os.path.join('data', 'fastas', 'probes', probe)):
                with open(os.path.join('data', 'fastas', 'probes', probe, file), 'r') as infile:
                    # Read the contents of the file and write them to the output file
                    outfile.write(infile.read())
                    print(f'Added {file} to output file.')
                    # Write a newline to ensure separation between files
                    outfile.write('\n')

    print(f"All FASTA entries have been concatenated into {output_file}")


def extractor_parser():
    Datahub.full_probes = table_parser()
    # download_genes(Datahub.full_probes)
    # concatenate_fastas()