# data/

```
Place the source input data within this directory path. Accepted formats are CSV/XLSX files and fields are listed below.
- (ensure that the file name begins with "source")
  * ex. ~/path/to/project/data/source_library.xlsx
  * ex. ~/path/to/project/data/source_lib.csv

<FIELDS>
- the default column names are listed in quotations.
- to make alterations to the column names, run <init --set_labels>. 

- required:
  > "SMILES"
    * (The compound in SMILES format)
  > "PGCC_score"
    * (PGCC post-treatment biodata)
    * [post_pgcc / pre_pgcc * 100]

- (optional)
  > "nonPGCC_score"
    * (nonPGCC post-treatment biodata)
    * [post_cancCell / pre_cancCell * 100]
  > "compound_name"
  > "pathway"
  > "target"
  > "info"
```