# Create hg38 annotation

Generates genome-wide gene and chromosome annotation datasets for the
human reference genome GRCh38 (hg38) using Ensembl.

This function retrieves all HGNC genes from Ensembl, annotates their
genomic positions and chromosome arm assignments using
[`get_genomic_info()`](get_genomic_info.md), and saves the resulting
datasets as RDS files.

This function was used to generate the packaged datasets hg38_gene_info
and hg38_chrom_info.

## Usage

``` r
create_hg38_annotation(
  gene_info_file = "hg38_gene_info.rds",
  chrom_info_file = "hg38_chrom_info.rds"
)
```

## Arguments

- gene_info_file:

  Output file for gene annotation.

- chrom_info_file:

  Output file for chromosome annotation.

## Value

A list containing the generated hg38_gene_info and hg38_chrom_info
files.
