# Variant Annotations Pipeline
Annotates a list of variants with:
* ClinVar clinical significance  
* gnomAD allele frequency  

## Setup

Create the environment:

```bash
conda env create -f environment.yml
conda activate <environment_name>
```

## Usage

```python main.py \
  --clinvar clinvar_ref_file/variant_summary.txt \
  --input input_files/variants.txt \
  --output results/output.tsv \
  --genome hg19
```


## Input Format

Plain text file with one variant per line:

```
1-94471075-A-G   
1-94473807-C-T
```

## Output

Tab-separated file containing:
* VariantID
* Clinical Significance 
* Review Status
* Number of Submissions
* Allele Frequencies

## Notes
* Genome Options: hg18 or hg38
* Uses gnomAD API (rate-limited, may run slowly)

# Variant_Annotations_ClinVar_gnomAD
This repo contains a script to run variant annotation on a file that contains list of variants.
