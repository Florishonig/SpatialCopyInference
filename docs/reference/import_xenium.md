# Import Xenium data

Imports Xenium spatial transcriptomics data using Seurat::LoadXenium()
and optionally adds cell type annotations from a CSV file.

## Usage

``` r
import_xenium(data_dir, annotation_file = NULL)
```

## Arguments

- data_dir:

  Path to Xenium directory containing default Xenium output files.

- annotation_file:

  Optional path to a CSV file containing cell type annotations.

## Value

A Seurat object containing Xenium data and associated metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
# Import Xenium data only
xenium.obj <- import_xenium(
  data_dir = "path/to/Xenium_outputs/outs",
  annotation_file = NULL)

# Import Xenium data and cell type annotations
xenium.obj <- import_xenium(
  data_dir = "path/to/Xenium_outputs/outs",
  annotation_file = "path/to/cell_type_annotations.csv")} # }
```
