import seq_utils
import utils
import defaults

from colored_logging import colored_logging
import logging

if __name__ == '__main__':
    colored_logging(log_file_name='ltr_fasta_blaster.txt')

    probe_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name='probe_dict.pkl')

    clean_tblastn2gb_results: dict = seq_utils.blast_retriever(object_dict=probe_dict,
                                                               command='tblastn',
                                                               genome=defaults.SPECIES,
                                                               online_database='nucleotide',
                                                               input_database_path=defaults.LTR_DB,
                                                               genbank_retrieval=True,
                                                               multi_threading=True)

    utils.pickler(data=clean_tblastn2gb_results,
                  output_directory_path=defaults.PICKLE_DIR,
                  output_file_name='ltr_fasta_blast.pkl')

    utils.directory_content_eraser(directory_path=defaults.TMP_DIR)

    logging.info(f'Successfully performed tBLASTn and retrieval on {len(clean_tblastn2gb_results)} sequences.')
