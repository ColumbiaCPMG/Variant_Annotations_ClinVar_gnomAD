#!/usr/bin/env python3

import argparse
import os


# ===================== CLI =============================

def parse_args():
    parser = argparse.ArgumentParser(description="Annotate variants with ClinVar + gnomAD AF")

    parser.add_argument("--input", required=True, help="Input file located in input_files/")
    parser.add_argument("--genome", required=True, help="hg19/GRCh37 or hg38/GRCh38")

    parser.add_argument("--sleep", type=float, default=7.0)
    parser.add_argument("--download-clinvar", action="store_true")

    return parser.parse_args()


# ===================== LIBRARIES =======================

import pandas as pd
import requests
import time
import gzip
import shutil


# ================= PATH SETUP =========================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

INPUT_DIR = os.path.join(BASE_DIR, "input_files")
CLINVAR_DIR = os.path.join(BASE_DIR, "clinvar_ref_file")
RESULTS_DIR = os.path.join(BASE_DIR, "results")


# ================= CLINVAR DOWNLOAD ====================

def download_clinvar_file(url):
    print("Downloading latest ClinVar variant summary...")

    os.makedirs(CLINVAR_DIR, exist_ok=True)

    gz_file = os.path.join(CLINVAR_DIR, "variant_summary.txt.gz")
    out_file = os.path.join(CLINVAR_DIR, "variant_summary.txt")

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(gz_file, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    print("Download complete. Unzipping...")

    with gzip.open(gz_file, "rb") as f_in:
        with open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    print(f"ClinVar file ready: {out_file}")
    return out_file


# ================= CLINVAR LOADING =====================

def load_clinvar(clinvar_file):

    print(f"Reading ClinVar file: {clinvar_file}")
    
    cols_to_use = [
        "Chromosome",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
        "ClinicalSignificance",
        "ReviewStatus",
        "NumberSubmitters"
    ]

    clinvar_df = pd.read_csv(
        clinvar_file,
        sep="\t",
        usecols=cols_to_use,
        dtype=str,
        low_memory=True,
        engine="c"
    )

    clinvar_df.columns = clinvar_df.columns.str.strip().str.replace("#", "")

    clinvar_df["VariantID"] = (
        clinvar_df["Chromosome"] + "-" +
        clinvar_df["PositionVCF"] + "-" +
        clinvar_df["ReferenceAlleleVCF"] + "-" +
        clinvar_df["AlternateAlleleVCF"]
    )

    return clinvar_df[
        ["VariantID", "ClinicalSignificance", "ReviewStatus", "NumberSubmitters"]
    ].dropna()


# ================= INPUT PROCESSING =====================

def load_input(input_file):
    file_path = os.path.join(INPUT_DIR, input_file)

    print(f"Reading input file: {file_path}")

    input_df = pd.read_csv(file_path, dtype=str)
    cols = [c.strip().lower() for c in input_df.columns]

    # detect header-like column
    for i, col in enumerate(cols):
        if any(k in col for k in ["variant", "var", "id"]):
            input_df = input_df.rename(columns={input_df.columns[i]: "VariantID"})
            break
    else:
        # no header case
        first_val = str(input_df.columns[0]).strip()

        if first_val and first_val[0].isdigit():
            print("No header detected. Treating file as variant list.")
            input_df = pd.read_csv(file_path, dtype=str, header=None)
            input_df.columns = ["VariantID"]
        else:
            input_df = input_df.rename(columns={input_df.columns[0]: "VariantID"})

    input_df["VariantID"] = (
        input_df["VariantID"]
        .astype(str)
        .str.replace('"', '', regex=False)
        .str.strip()
    )

    return input_df[input_df["VariantID"] != ""]


# ================= MERGE ===============================

def merge_data(input_df, clinvar_subset):
    return input_df.merge(clinvar_subset, on="VariantID", how="left")


# ================= DATASET =============================

def select_dataset(genome_build):
    gb = genome_build.lower().strip()

    if gb in ["hg19", "grch37"]:
        return "gnomad_r2_1"
    elif gb in ["hg38", "grch38"]:
        return "gnomad_r3"
    else:
        raise ValueError("Genome must be hg19/GRCh37 or hg38/GRCh38")


# ================= API ================================

def fetch_block(variant_id, dataset, block, retries=3):
    url = "https://gnomad.broadinstitute.org/api"

    query = {
        "query": f"""
        {{
          variant(variantId: "{variant_id}", dataset: {dataset}) {{
            {block} {{
              af
              ac
              an
            }}
          }}
        }}
        """
    }

    for _ in range(retries):
        try:
            r = requests.post(url, json=query)
            data = r.json()

            if "errors" in data:
                return None

            variant = data.get("data", {}).get("variant")
            return variant.get(block) if variant else None

        except:
            time.sleep(0.5)

    return None


def fetch_af(variant_id, dataset):
    genome = fetch_block(variant_id, dataset, "genome")
    exome = fetch_block(variant_id, dataset, "exome")

    genome_af = genome.get("af") if genome else None
    exome_af = exome.get("af") if exome else None

    genome_ac = genome.get("ac") if genome else 0
    genome_an = genome.get("an") if genome else 0

    exome_ac = exome.get("ac") if exome else 0
    exome_an = exome.get("an") if exome else 0

    total_ac = (genome_ac or 0) + (exome_ac or 0)
    total_an = (genome_an or 0) + (exome_an or 0)

    total_af = total_ac / total_an if total_an else None

    return {
        "genome_af": genome_af,
        "exome_af": exome_af,
        "total_af": total_af
    }


# ================= ANNOTATION ==========================

def annotate_af(annotation_df, dataset, sleep_time):
    results = []
    total = len(annotation_df)

    for i, vid in enumerate(annotation_df["VariantID"], start=1):
        results.append(fetch_af(vid, dataset))

        if i % 10 == 0 or i == total:
            print(f"Processed {i}/{total} variants")

        time.sleep(sleep_time)

    af_df = pd.DataFrame(results)

    return pd.concat(
        [annotation_df.reset_index(drop=True), af_df],
        axis=1
    )


# ================= OUTPUT ==============================

def save_output(df, input_file):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(RESULTS_DIR, f"{base_name}_results.csv")

    df.to_csv(output_file, index=False)
    print(f"Saved results to: {output_file}")


# ================= PIPELINE ============================

def run_pipeline(input_file, genome_build, sleep_time, download_flag):
    print("\n=== Starting Variant Annotation Pipeline ===\n")

    clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    clinvar_file = os.path.join(CLINVAR_DIR, "variant_summary.txt")

    if download_flag:
        print("[1/5] Downloading ClinVar reference...")
        clinvar_file = download_clinvar_file(clinvar_url)
    else:
        print(f"[1/5] Using existing ClinVar file: {clinvar_file}")

    print("[2/5] Loading ClinVar data...")
    clinvar_subset = load_clinvar(clinvar_file)
    print(f"Loaded {len(clinvar_subset):,} ClinVar records")

    print("[3/5] Processing input file...")
    input_df = load_input(input_file)
    print(f"Loaded {len(input_df):,} variants")

    print("[4/5] Merging with ClinVar...")
    annotation_df = merge_data(input_df, clinvar_subset)

    dataset = select_dataset(genome_build)
    print(f"[5/5] Fetching gnomAD data using dataset: {dataset}")
    print("This may take a while depending on number of variants.\n")

    annotation_df = annotate_af(annotation_df, dataset, sleep_time)

    save_output(annotation_df, input_file)

    print("\n=== Pipeline Complete ===\n")


# ================= MAIN ================================

def main():
    args = parse_args()

    run_pipeline(
        input_file=args.input,
        genome_build=args.genome,
        sleep_time=args.sleep,
        download_flag=args.download_clinvar
    )


if __name__ == "__main__":
    main()