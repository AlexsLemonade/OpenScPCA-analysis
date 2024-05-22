#!/usr/bin/env python3

import anndata
import scrublet as scr
from pathlib import Path
import pyprojroot

# Function to run scrublet
def run_scrublet(adata: anndata) -> tuple[float, str]:
    """
        Run scrublet on a counts matrix from an AnnData object, following https://github.com/swolock/scrublet?tab=readme-ov-file#quick-start
        Returns a tuple of (barcodes, doublet_scores, predicted_doublets)
    """

    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # stringify
    doublet_scores_str = [str(x) for x in doublet_scores]
    predicted_doublets_str = [("doublet" if x else "singlet") for x in predicted_doublets]

    return( (list(adata.obs_names), doublet_scores_str, predicted_doublets_str) )

def main() -> None:

    # Define directories
    openscpca_base = pyprojroot.find_root(pyprojroot.has_dir(".git"))
    module_base = Path(openscpca_base / "analyses/doublet-detection")

    data_dir = Path(module_base / "scratch/benchmark_datasets")
    result_dir = Path(module_base / "results/benchmark_results")
    result_dir.mkdir(parents = True, exist_ok = True)

    datanames = ["hm-6k", "pbmc-1B-dm", "pdx-MULTI", "HMEC-orig-MULTI"]

    for dataname in datanames:
        input_anndata = dataname + "_anndata.h5ad"
        result_tsv = dataname + "_scrublet.tsv"

        adata = anndata.read_h5ad( Path(data_dir / input_anndata) )
        scrub_results = run_scrublet(adata)

        with open(Path(result_dir / result_tsv), "w") as f:
            f.write("barcode\tdoublet_score\tpredicted_doublets\n")
            for x in range(len(scrub_results[0])):
                f.write("\t".join( [scrub_results[0][x], scrub_results[1][x], scrub_results[2][x]] ) + "\n")

main()