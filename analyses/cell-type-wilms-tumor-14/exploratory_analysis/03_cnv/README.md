## Exploratory analysis

#### CNV analysis
The goal of this exploratory CNV analysis is to see if `CopyKat` and `inferCNV` can be useful methods for predicting tumor cells based on ploidy. This analysis uses the library `SCPCL000850` (sample `SCPCS000518`) to explore these methods.

See below for details on how to run this exploratory CNV analysis.
```bash
cd /path/to/OpenScPCA-analysis
cd analyses/cell-type-wilms-tumor-14/exploratory_analysis/03_cnv
# required for rendering Rmd files
export RSTUDIO_PANDOC="/usr/lib/rstudio/resources/app/bin/quarto/bin/tools"

## CopyKat
# generate CNV results using copykat
Rscript ./03_runCopyKat.R
# render Rmd
Rscript -e "rmarkdown::render('./03_copykat_SCPCL000850.Rmd')"

## inferCNV
# generate inferCNV inputs
step_name="03_cnv"
scratch_dir_step="../../scratch/${step_name}/infercnv" && mkdir -p ${scratch_dir_step}
Rscript ../../../cell-type-ewings/scripts/cnv-workflow/00-make-gene-order-file.R \
  --local_ref_dir ${scratch_dir_step} \
  --scratch_dir ${scratch_dir_step}
# run inferCNV (takes a while)
Rscript ./03_runInferCNV.R
# render Rmd
Rscript -e "rmarkdown::render('./03_infercnv_SCPCL000850.Rmd')"
```