## Exploratory analysis
```bash
cd /path/to/OpenScPCA-analysis
cd analyses/cell-type-wilms-tumor-14
# required for rendering Rmd files
export RSTUDIO_PANDOC="/usr/lib/rstudio/resources/app/bin/quarto/bin/tools"

## CopyKat
# generate CNV results using copykat
Rscript exploratory_analysis/03_runCopyKat.R
# render Rmd
Rscript -e "rmarkdown::render('exploratory_analysis/03_copykat_SCPCL000850.Rmd')"

## inferCNV
# generate inferCNV inputs
step_name="03_cnv"
scratch_dir_step="scratch/${step_name}/infercnv" && mkdir -p ${scratch_dir_step}
Rscript ../cell-type-ewings/scripts/cnv-workflow/00-make-gene-order-file.R \
  --local_ref_dir ${scratch_dir_step} \
  --scratch_dir ${scratch_dir_step}
# run inferCNV (takes a while)
Rscript exploratory_analysis/03_runInferCNV.R
# render Rmd
Rscript -e "rmarkdown::render('exploratory_analysis/03_infercnv_SCPCL000850.Rmd')"
```