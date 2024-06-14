import pandas as pd
import numpy as np
import os
from collections import defaultdict
from Bio import Entrez, SeqIO
import re

def seq_regex_parser():
    # Define the regular expression pattern
    pattern = r'(?P<accession_number>[A-Z]{2}_\d{6}\.\d)'

    # Find all matches in the sequence
    matches = re.finditer(pattern, sequence)

    # Create a dictionary to store the results
    results = defaultdict(list)

    # Iterate over the matches and extract the information
    for match in matches:
        accession_number = match.group('accession_number')
        results['accession_numbers'].append(accession_number)

    # Retrieve the last 200 characters of the sequence
    last_200 = sequence[-200:]
    results['last_200'] = last_200


    return results