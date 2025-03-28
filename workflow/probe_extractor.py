from colored_logging import colored_logging
import pandas as pd
import logging
import os

from object_class import Object
import seq_utils
import utils
import defaults


def table_parser(input_csv_file):
    """
    Parse the CSV files containing the probe genes and their respective IDs

        Parameters
        ----------
            :param input_csv_file: The path to the CSV file containing probe info.

        Returns
        -------
            :returns: A dictionary containing object pairs.
    """
    # Read the CSV file
    probe_table = pd.read_csv(input_csv_file)

    probe_dict: dict = {
        f'{str(row["Accession"])}': Object(
            family=str(row['Family']),
            virus=str(row['Name']),
            abbreviation=str(row['Abbreviation']),
            probe=str(row['Probe']),
            accession=str(row['Accession']),
            identifier=utils.random_string_generator(6)
        )
        for index, row in probe_table.iterrows()
    }
    return probe_dict


if __name__ == '__main__':
    colored_logging(log_file_name='probe_extractor.txt')

    probe_dict: dict = table_parser(input_csv_file=os.path.join(defaults.TABLE_INPUT_DIR, 'probes.csv'))

    probe_extraction: dict = seq_utils.gb_executor(object_dict=probe_dict,
                                                   online_database='protein',
                                                   display_full_info=defaults.DISPLAY_OPERATION_INFO,
                                                   display_warning=defaults.DISPLAY_REQUESTS_WARNING)

    utils.pickler(data=probe_extraction,
                  output_directory_path=defaults.PICKLE_DIR,
                  output_file_name='probe_dict.pkl')