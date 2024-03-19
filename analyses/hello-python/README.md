# Hello Python

This is an example analysis module using a Jupyter Notebook to perform a simple analysis of the OpenScPCA data.

The analysis counts the number of cells in each processed SCE file in the current data directory.

## Usage

To run the analysis, you will first need to activate the included conda environment.
From this analysis directory, run the following commands at the command line:

```bash
conda env update -f environment.yml -n openscpca-hello-python
conda activate openscpca-hello-python
```

Then run the following command from this analysis directory:

```bash
bash hello-python.sh
```

### Output

- `hello.html`: a rendered version of the Jupyter Notebook
- `results/cell_counts.csv`: a table of counts for each library
- `plots/cell_counts.pdf`: histograms of cell counts by project
