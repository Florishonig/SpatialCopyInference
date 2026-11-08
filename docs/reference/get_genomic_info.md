# Get genomic information

Retrieves genomic coordinates and chromosome arm annotations for a set
of genes from Ensembl and combines them with chromosome size and
centromere information.

The function returns gene level and chromosome level annotation tables
that can be used for chromosome arm spatial copy number inference.

This function was used to generate the packaged datasets hg38_gene_info
and hg38_chrom_info.

## Usage

``` r
get_genomic_info(
  xenium.obj = NULL,
  genes = NULL,
  chrom_sizes_file = system.file("extdata", "hg38.chrom.sizes.txt", package =
    "SpatialCopyInference"),
  centromere_file = system.file("extdata", "centromeres_hg38.tsv", package =
    "SpatialCopyInference")
)
```

## Arguments

- xenium.obj:

  Xenium object containing genes to annotate.

- genes:

  Character vector of HGNC gene symbols. Either `xenium.obj` or `genes`
  must be supplied.

- chrom_sizes_file:

  Path to chromosome size file.

- centromere_file:

  Path to centromere annotation file.

## Value

A list containing gene and chromosome annotation information.
