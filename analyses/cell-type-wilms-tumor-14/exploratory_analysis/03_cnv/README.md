## Exploratory analysis

#### CNV analysis
The goal of this exploratory CNV analysis is to see if `CopyKat` and `inferCNV` can be useful methods for predicting tumor cells based on ploidy. This analysis uses the library `SCPCL000850` (sample `SCPCS000518`) to explore these methods.

See below for details on how to run this exploratory CNV analysis. 
- Note: Outputs from main pipeline step `00_preprocessing_rds` for library `SCPCL000850` are required to run the following analysis, see [scripts](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/analyses/cell-type-wilms-tumor-14/run_cell-type-wilms-14.sh) to run main pipeline.
```bash
# run this analysis from analysis module folder to make use of renv.lock
cd /path/to/OpenScPCA-analysis
cd analyses/cell-type-wilms-tumor-14/

# required for rendering Rmd files.
# This path exists on OpenScPCA Lightsail for Research instances, but you may need to specify a different path if running on a different platform, or run `rmarkdown::render` inside your rstudio session.
export RSTUDIO_PANDOC="/usr/lib/rstudio/resources/app/bin/quarto/bin/tools"

## CopyKat
# generate CNV results using copykat 
Rscript ./exploratory_analysis/03_cnv/03_runCopyKat.R
# render Rmd
Rscript -e "rmarkdown::render('./exploratory_analysis/03_cnv/03_copykat_SCPCL000850.Rmd')"

## inferCNV
# generate inferCNV inputs
step_name="03_cnv"
scratch_dir_step="./scratch/${step_name}/infercnv" && mkdir -p ${scratch_dir_step}
Rscript ../cell-type-ewings/scripts/cnv-workflow/00-make-gene-order-file.R \
  --local_ref_dir ${scratch_dir_step} \
  --scratch_dir ${scratch_dir_step}
# run inferCNV (takes a while)
Rscript ./exploratory_analysis/03_cnv/03_runInferCNV.R
# render Rmd
Rscript -e "rmarkdown::render('./exploratory_analysis/03_cnv/03_infercnv_SCPCL000850.Rmd')"
```