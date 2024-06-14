from extractor_parser import extractor_parser, table_parser, download_genes, concatenate_fastas

if __name__ == '__main__':
    gene_dicts = table_parser()
    # download_genes(gene_dicts)
    # concatenate_fastas()