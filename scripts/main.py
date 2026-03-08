import os
import argparse
import subprocess
import glob
import pandas as pd
from tqdm import tqdm
from itertools import product
from concurrent.futures import ThreadPoolExecutor, as_completed


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run PlasmidNL typing pipeline (in parallel)"
    )

    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input directory containing plasmid FASTA files",
    )

    parser.add_argument(
        "--identity",
        "-id",
        type=int,
        default=90,
        help="Identity value cut-off for detecting genes (default: 90)",
    )

    parser.add_argument(
        "--coverage",
        "-cv",
        type=int,
        default=80,
        help="Identity value cut-off for detecting genes (default: 80)",
    )

    parser.add_argument(
        "--jobs", "-n", type=int, default=1, help="Number of parallel jobs (default: 1)"
    )

    parser.add_argument(
        "--keep_intermediate_files",
        "-ki",
        action="store_true",
        help="Keep intermediate files",
    )

    return parser.parse_args()


def run_blast_pair(pair: tuple, identity: int, in_dir: str, target_dir: str) -> str:
    """
    Function to run a single blast between a pair of FASTAs

    Parameters
    ----------
    pair : tuple
        tuple containing each member of the pair as an element
    identity : int
        threshold to apply on the hit identity value
    in_dir : str
        directory where the first FASTA file is stored
    target_dir : str
        directory where the second FASTA file is stored

    Returns
    -------
    str
        output message for which pair was run
    """
    # Seperate the pair
    query, subject = pair
    query_base = os.path.splitext(os.path.basename(query))[0]
    subject_base = os.path.splitext(os.path.basename(subject))[0]

    # Prepare directories
    intermediate_dir = os.path.join("output", "intermediate", query_base)
    os.makedirs(intermediate_dir, exist_ok=True)
    log_dir = os.path.join("logs", query_base)
    os.makedirs(log_dir, exist_ok=True)

    # File paths
    out_file = os.path.join(intermediate_dir, f"{subject_base}.csv")
    log_path = os.path.join(log_dir, f"{subject_base}.log")

    # BLAST command
    outfmt = "10 qseqid sseqid qlen slen qstart qend sstart send evalue score length pident mismatch gapopen"
    cmd = [
        "blastn",
        "-query",
        os.path.join(in_dir, query),
        "-subject",
        os.path.join(target_dir, subject),
        "-outfmt",
        outfmt,
        "-perc_identity",
        str(identity),
    ]

    # Run BLAST
    with open(out_file, "a") as out_f, open(log_path, "w") as log_f:
        subprocess.run(cmd, stdout=out_f, stderr=log_f, check=True)

    return f"Finished: {query_base} vs {subject_base}"


def merge_blast_outputs(max_workers: int = 4):
    """
    Merges all the individual blast output files into one

    Parameters
    ----------
    max_workers : int, optional
        number of parallel jobs, by default 4

    """
    output_file = "output/full_blast_output.csv"
    log_path = "logs/file_join.log"
    os.makedirs("logs", exist_ok=True)

    print("\nMerging output files in parallel\n")

    # Remove existing output file if it exists
    if os.path.exists(output_file):
        os.remove(output_file)

    # Header line
    header = "file_name,qseqid,sseqid,qlen,slen,qstart,qend,sstart,send,evalue,score,length,pident,mismatch,gapopen\n"
    with open(output_file, "w") as out_f:
        out_f.write(header)

    # Find all CSV files in subdirectories
    csv_files = glob.glob("output/intermediate/*/*.csv")

    # Sort by query directory name to preserve consistent output order
    csv_files.sort(key=lambda f: os.path.basename(os.path.dirname(f)))

    # Function to read a CSV file and prepend the dirname
    def read_csv_with_dir(file_path):
        dirname = os.path.basename(os.path.dirname(file_path))
        lines = []
        try:
            with open(file_path, "r") as f:
                for line in f:
                    lines.append(f"{dirname},{line}")
        except Exception as e:
            with open(log_path, "a") as log_f:
                log_f.write(f"Error processing {file_path}: {e}\n")
        return lines

    # Parallel reading
    all_lines = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(read_csv_with_dir, f): f for f in csv_files}
        for future in tqdm(
            as_completed(futures), total=len(futures), desc="Merging CSVs"
        ):
            all_lines.extend(future.result())

    # Write all lines to the final output
    with open(output_file, "a") as out_f:
        out_f.writelines(all_lines)

    print(f"\nMerged {len(csv_files)} files into {output_file}")


# Function to filter overlapping hits for each qseqid
def filter_overlaps(group: pd.DataFrame) -> pd.DataFrame:
    """
    Filters

    Parameters
    ----------
    group : pd.DataFrame
        dataframe containing all hits from one contig

    Returns
    -------
    pd.DataFrame
        dataframe filtered for overlapping hits
    """
    kept_hits = []
    for _, row in group.iterrows():
        overlap = False
        for kept in kept_hits:
            # Check if this hit overlaps with any kept hit
            if not (row["qend"] < kept["qstart"] or row["qstart"] > kept["qend"]):
                overlap = True
                break
        if not overlap:
            kept_hits.append(row)
    return pd.DataFrame(kept_hits)


def remove_directory_tree(start_directory: str):
    """Recursively and permanently removes the specified directory, all of its
    subdirectories, and every file contained in any of those folders."""
    for name in os.listdir(start_directory):
        path = os.path.join(start_directory, name)
        if os.path.isfile(path):
            os.remove(path)
        else:
            remove_directory_tree(path)
    os.rmdir(start_directory)


def maximum_variant(blast_df: pd.DataFrame) -> pd.DataFrame:
    # Ensure we don't modify the original DataFrame
    df = blast_df.copy()

    # Sort to ensure tie-breaking order: lowest evalue, highest pident
    df_sorted = df.sort_values(
        by=["Contig", "evalue", "pident"], ascending=[True, True, False]
    )

    # Keep only the top match for each Contig and sseqid
    top_matches = df_sorted.groupby(["Contig", "sseqid"], as_index=False).first()

    # For each Contig-sseqid, get the max pident (vectorized)
    max_pident = top_matches.groupby(["Contig", "sseqid"], as_index=False)[
        "pident"
    ].max()

    # Pivot to get sseqid as columns, Contig as rows
    gene_summary = max_pident.pivot(
        index="Contig", columns="sseqid", values="pident"
    ).reset_index()

    # Save to CSV
    gene_summary.to_csv("output/variant_identity_table.csv", sep=";", index=False)

    # Copy and prepare the DataFrame
    gene_summary_index = gene_summary.set_index("Contig")

    # Identify gene groups from column names
    # genes = [col for col in blast_output.columns]
    genes = os.listdir("targets")
    genes = [x.split(".fa")[0] for x in genes]
    gene_groups = list({col.split("_")[0].split("(")[0] for col in genes})

    existing_cols = set(gene_summary_index.columns)

    group_dfs = []

    for gene_group in gene_groups:

        # genes belonging to this family
        group_cols = [
            col for col in genes if col.split("_")[0].split("(")[0] == gene_group
        ]

        # keep only columns that exist in dataframe
        present_cols = [c for c in group_cols if c in existing_cols]

        if present_cols:
            group_df = gene_summary_index[present_cols]

            max_vals = group_df.max(axis=1)

            valid_rows = group_df.notna().any(axis=1)
            max_cols = pd.Series(index=group_df.index, dtype=object)
            max_cols.loc[valid_rows] = group_df.loc[valid_rows].idxmax(axis=1)
        else:
            # gene family not found in any sequence
            max_vals = pd.Series(index=gene_summary_index.index, dtype=float)
            max_cols = pd.Series(index=gene_summary_index.index, dtype=object)

        summary_df = pd.DataFrame(
            {f"{gene_group}": max_cols, f"{gene_group}_identity": max_vals},
            index=gene_summary_index.index,
        )

        group_dfs.append(summary_df)

    all_genes = pd.concat(group_dfs, axis=1)
    all_genes = all_genes.reindex(sorted(all_genes.columns), axis=1)

    all_genes.to_csv("output/maximum_variant.csv", sep=";")

    # ---- presence/absence matrix ----
    # identity_cols = [c for c in all_genes.columns if c.endswith("_identity")]
    gene_cols = [c for c in all_genes.columns if not c.endswith("_identity")]

    all_genes_names = all_genes[gene_cols]

    out_df = all_genes_names.notna().astype(int)
    out_df.to_csv("output/presence_absence_table.csv")


def blast_files():
    """
    Main function that runs the blast and filtering.
    """
    # Collect input arguments
    args = parse_args()
    in_dir = args.input
    id_thresh = args.identity
    cov_thresh = args.coverage
    n_jobs = args.jobs
    keep_files = args.keep_intermediate_files

    # Possible file extensions that can be run
    approved_extensions = (".fasta", ".fa", ".fas", ".fna", ".ffn", ".faa")

    # Collect queries and targets
    queries = [x for x in os.listdir(in_dir) if x.endswith(approved_extensions)]
    targets = [x for x in os.listdir("targets") if x.endswith(approved_extensions)]

    all_pairs = list(product(queries, targets))
    job_total = len(all_pairs)

    print("\n\tChecking gene presence in input FASTAs...\n")

    # Run BLAST in parallel
    results = []
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        future_to_pair = {
            executor.submit(run_blast_pair, pair, id_thresh, in_dir, "targets"): pair
            for pair in all_pairs
        }
        for future in tqdm(as_completed(future_to_pair), total=job_total):
            pair = future_to_pair[future]
            try:
                results.append(future.result())
            except Exception as e:
                print(f"Error running BLAST for {pair}: {e}")

    print(f"\nFinished {job_total} BLAST jobs.")

    # Merge all the outputs
    merge_blast_outputs(max_workers=args.jobs)

    # Filter out the bad hits
    print("\tRemoving spurious hits\n")
    blast_out = pd.read_csv("output/full_blast_output.csv", sep=",")
    blast_thres = blast_out.loc[blast_out["pident"] >= id_thresh]
    blast_thres["scov"] = (
        (abs(blast_thres["sstart"] - blast_thres["send"]) + 1) / blast_thres["slen"]
    ) * 100
    blast_thres = blast_thres.loc[blast_thres["scov"] >= cov_thresh]
    blast_thres = blast_thres.sort_values(
        by=["qseqid", "pident"], ascending=[True, False]
    )

    # Create a top blast hit file
    blast_filtered = blast_thres.groupby("qseqid", group_keys=False).apply(
        filter_overlaps, include_groups=False
    )

    # Make column names more human readable
    blast_filtered = blast_filtered.rename(
        columns={
            "file_name": "Fasta",
            "qseqid": "Contig",
            "sseqid": "Gene",
            "qlen": "Hit length on contig",
            "slen": "Hit length on gene",
            "qstart": "Start position on contig",
            "qend": "End position on contig",
            "sstart": "Start position on gene",
            "send": "End position on contig",
            "length": "Hit length total",
            "pident": "Hit identity (%)",
            "mismatch": "Mismatched basepairs",
            "gapopen": "Hit gaps",
            "scov": "Gene coverage",
        }
    )

    # Save and remove intermediate files
    blast_filtered.to_csv("output/top_hits_raw.csv", index=False)
    if not keep_files:
        print("Removing intermediate files\n")
        remove_directory_tree("output/intermediate")

    # Create maximum variant file
    maximum_variant(blast_thres.rename(columns={"qseqid": "Contig"}))

    print(
        """
        Finished annotating genes associated with hypervirulance in Klebsiella pneumoniae

        Thank you for using this pipeline!

        If you use this pipeline in your research, please cite the following publication:

        Vendrik KEW, Teunis G, Anema F, Landman F, de Haan A, Bos J, Witteveen S,
        Schoffelen AF, de Greeff SC, Kuijper EJ, Hendrickx APA, Notermans DW;
        hvKp study group. Genomic epidemiology of putative hypervirulent Klebsiella
        pneumoniae species complex in Dutch patients, January-December 2022.
        Microbiol Spectr. 2026 Feb 3;14(2):e0225925. doi: 10.1128/spectrum.02259-25.
        """
    )


if __name__ == "__main__":
    blast_files()
