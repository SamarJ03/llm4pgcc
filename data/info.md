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

```
PATH LAYOUT PLAN FOR DATA/:

source_lib.xlsx (user input csv or xlsx file)
raw_lib.csv (user input in CSV format, only required params kept)
clean_lib.csv (raw_lib.csv after data preprocessing steps)
lib.json (JSON file of clean_lib.csv compounds with BOTH required and optional fields)
features.csv (CSV of all compounds from clean_lib.csv and maccs/ecfp4/rdkit descriptors)
features/
 - ecfp4_meta.csv (ecfp4 descriptors prior to low-var removal)
 - maccs_meta.csv (maccs descriptors prior to low-var removal)
 - rdkit_meta.csv (rdkit descriptors prior to low-var removal)
 - ridLowVar/ (low var descriptors removed)
   - ecfp4.csv
   - maccs.csv
   - rdkit.csv
scaffold/
 - ecfp4/ (contains train/test/valid sets)
 - maccs/ (contains train/test/valid sets)
 - rdkit/ (contains subgroupings train/test/valid sets including file of all rdkit descriptors train/test/valid set)
```