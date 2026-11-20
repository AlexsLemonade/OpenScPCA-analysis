* `01-convert-osteocar.R` converts a given input Seurat reference object to SCE and AnnData objects.

* `02-create-reference.R` takes the primary and metastasis OsteoCAR references (SCE objects) and creates a single harmonized reference with one level of annotation where Tumor and Mesenchymal labels come from Ann_Level1, and all other (normal) cell type labels come from Ann_Level2. 
After labels are harmonized, each reference is aggregated independently and then combined to create a single reference.

* `03-train-singler-model.R` trains a SingleR model from the harmonized, aggregated reference.
An ScPCA SCE object is used to define the genes the model expects in query data, after converting its Ensembl IDs to gene symbols with `rOpenScPCA::sce_to_symbols()`.
Query SCE objects must be converted the same way before classification.
The trained model is saved to `scratch/`.