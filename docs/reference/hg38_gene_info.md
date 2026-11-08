# hg38 gene information

Gene level information for the human reference genome GRCh38 (hg38).

## Usage

``` r
hg38_gene_info
```

## Format

A data frame with 41,480 rows and 24 columns containing Ensembl gene
identifiers, HGNC symbols, genomic coordinates, chromosome assignments
and chromosome arm annotations.

## Source

Ensembl

## Details

The dataset was generated from Ensembl gene annotations using the
SpatialCopyInference functions
[`get_genomic_info()`](get_genomic_info.md) and
[`create_hg38_annotation()`](create_hg38_annotation.md).

It contains genomic coordinates, chromosome assignments and chromosome
arm annotations for all HGNC genes annotated in Ensembl.
