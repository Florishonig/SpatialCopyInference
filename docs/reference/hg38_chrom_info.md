# hg38 chromosome information

Chromosome and chromosome arm information for the human reference genome
GRCh38 (hg38).

## Usage

``` r
hg38_chrom_info
```

## Format

A data frame with 23 rows and 13 columns containing chromosome lengths,
centromere positions, chromosome arm boundaries and cumulative genomic
coordinates.

## Source

Ensembl

## Details

The dataset was generated from Ensembl chromosome and centromere
coordinates using the SpatialCopyInference functions
[`get_genomic_info()`](get_genomic_info.md) and
[`create_hg38_annotation()`](create_hg38_annotation.md).
