* `01-convert-osteocar.R` converts a given input Seurat reference object to SCE and AnnData objects.

* `02-create-reference.R` takes the primary and metastasis OsteoCAR references (SCE objects) and creates a single harmonized reference with one level of annotation where Tumor and Mesenchymal labels come from Ann_Level1, and all other (normal) cell type labels come from Ann_Level2. 
After labels are harmonized, each reference is aggregated independently and then combined to create a single reference.

* `03-train-singler-model.R` trains a SingleR model from the harmonized, aggregated reference.
An ScPCA SCE object is used to define the genes the model expects in query data, after converting its Ensembl IDs to gene symbols with `rOpenScPCA::sce_to_symbols()`.
Query SCE objects must be converted the same way before classification.
The trained model is saved to `scratch/`.

* `04-annotate-sample-singler.R` annotates a single ScPCA SCE object using the trained SingleR model from `03-train-singler-model.R`.
The SCE's Ensembl IDs are converted to gene symbols with `rOpenScPCA::sce_to_symbols(reference = "10x2020")`, the same conversion used when training the model.
The script saves the full SingleR result as an RDS file and a TSV file with the cell barcodes, labels, and `delta.next` values.
Input and output paths are provided as arguments, so the script can be run on each sample separately.