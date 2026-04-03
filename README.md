# Variant Annotation Pipeline
 
The pipeline annotates a list of genetic variants with:
* ClinVar clinical significance  
* gnomAD allele frequency  
 
## Setup

Clone the repository locally. The project primarily relies on standard Python libraries along with a few common external dependencies.
If you encounter missing modules or installation issues, you can create a dedicated environment using the provided environment.yml file:
 
```
conda env create -f environment.yml
conda activate <environment_name>
```
## ClinVar Reference File
The ClinVar variant summary file can be downloaded from the official NCBI FTP source:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
Alternatively, the pipeline can automatically download this file using the --download-clinvar flag (see Usage section).
 
## Usage

From the main project folder use the following command to run the pipeline:
```
python3 scripts/annotate_variants.py --input <input_file> --genome hg19
```
Required arguments:  
--input: Name of the input file located in the input_files/ directory  
--genome: Reference genome build (hg19/GRCh37 or hg38/GRCh38)

Optional arguments:  
--download-clinvar: Downloads the latest ClinVar variant summary file before running the pipeline --sleep: Time delay (in seconds) between API requests (default: 7.0)

* By default the results file will have the same name as the input file with “_results” suffix at the end in the results folder.
* The sleep timer between API queries is 7.0 seconds which is suggested as a smaller timer might send too many requested and get the IP blocked.
* Also, use —download-clinvar the first time you run it, to download the latest clinvar file. After initial run only use it when you want to update the clinvar files for subsequent analysis.
* Output files are automatically saved in the results/ directory.
The output filename is derived from the input filename with a _results.csv suffix
Example:
```
test_variants.txt → test_variants_results.csv
```
 
## Input Format
Plain text file with one variant per line:
```
1-94471075-A-G   
1-94473807-C-T
```
 
## Output
The pipeline generates a tabular file containing the following fields:  
* VariantID
* Clinical Significance
* Review Status
* Number of Submissions
* Allele Frequencies (gnomAD)

## Notes
Supported genome builds:  
* hg19 / GRCh37
* hg38 / GRCh38

The pipeline uses the gnomAD API, which is rate-limited. Execution may take time depending on the number of variants.  

In some cases, the gnomAD API may not return results for certain variants even if they exist in the browser. This is often due to differences in:  
* genome build
* variant representation (normalization)
* dataset availability

It is recommended to manually verify such variants on the gnomAD browser if needed.

## Date
April 3, 2026

## Author
Mallika
