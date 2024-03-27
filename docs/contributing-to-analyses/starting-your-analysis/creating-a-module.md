# Creating an analysis module

Once you have [discussed your proposed analysis](../../communications-tools/index.md#github-discussions) and [filed an issue to get started](../../communications-tools/index.md#github-issues), you're ready to create an analysis module.

We have provided a script [`create-analysis-module.py`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/create-analysis-module.py) you can use to make this process easier.
Before you begin, please review [this documentation](../analysis-modules/index.md) for an overview of skeleton file structure of an analysis module.


!!! note "Is this your first analysis module?"
    If you are creating your first analysis module, we recommend that you [scope your first pull request](../creating-pull-requests/scoping-pull-requests.md)  to contain only:

    - The analysis module skeleton created by running `create-analysis-module.py`
    - Some [initial documentation of your module](STUB_LINK documenting analysis modules) in the `README.md` file

    This way, you can get some experience filing pull requests, undergoing code review, and using other Git skills before you really get going on the analysis.


## Running the module creation script

You should use the [`create-analysis-module.py`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/create-analysis-module.py) script to create a skeleton module for your analysis in the `analyses` folder.

Before running this script, you should decide whether you will use R or Python (or both!) in your module.
When running the script, you can use flags to add additional files and/or folders to your module skeleton that you will need for your chosen language, [as described below](#script-flags).


To run this script, take these steps:

1. Open a [terminal](../../software-platforms/general-tools/using-the-terminal.md) window.
    - You may wish to launch terminal [from GitKraken](../../software-platforms/general-tools/using-the-terminal.md#gitkraken) so that you are automatically placed in the repository.

1. Make sure you are working in your `base` conda environment by running `conda activate base`.

1. As needed, navigate using `cd` to your `OpenScPCA-analysis` fork.

1. Run the module creation script using your desired flags:

```bash
# Create a module called `my-module-name`, for example
./create-analysis-module.py my-module-name <additional flags here>
```

!!! tip
    You can always run the script with the `--help` flag to reveal the help menu and see other flags you can use:

    ```bash
    ./create-analysis-module.py --help
    ```

## Module creation script flags


We recommend using one of these flags when creating your module.
Each flag will create a skeleton module with the given files and folders.

Note that you can use both an R and a Python flag if you want to write your module in both languages.

### Flags to create an R module

#### `--use-r`

This flag will add a template R notebook to your module.
<div class="grid" markdown>

You can use this flag as:

```bash
# Create a module called `my-module-name` with an example R Markdown notebook
./create-analysis-module.py my-module-name --use-r
```

This will create skeleton module called `my-module-name` with the files and folders shown.

```{ .no-copy title="Module directory with --use-r flag"}
├── scripts
│   └── ...
├── results
│   └── README.md
├── plots
│   └── ...
├── scratchIG
│   └── ...
├── README.md
├── notebook-template.Rmd
├── Dockerfile
├── .gitignore
└── .dockerignore
```
</div>

- You can use `notebook-template.Rmd` as a starting point for any R Markdown notebooks you create while writing your analysis



#### `--use-renv`

!!! info
    The `hello-r` example module was created with this flag.
    Learn more about [using `renv` to manage your R environment](../starting-your-analysis/determining-software-requirements/#using-renv).

- This flag will both:
    - Add a template R notebook to your module
    - Initialize an `renv` environment for your module

You can use this flag as:

```bash
# Create a module called `my-module-name` with an example R Markdown notebook and an `renv` environment
./create-analysis-module.py my-module-name --use-renv
```

This will create skeleton module called `my-module-name` with these files and folders:


```
├── scripts
│   └── ...
├── results
│   └── README.md
├── plots
│   └── ...
├── scratch
│   └── ...
├── README.md
├── notebook-template.Rmd
├── components
│   └── dependencies.R
├── renv.lock
├── renv
│   └── ...
├── .Rprofile
├── Dockerfile
├── .gitignore
└── .dockerignore
```

- You can use `notebook-template.Rmd` as a starting point for any R Markdown notebooks you create while writing your analysis
- You can use `components/dependencies.R` to [pin R package dependencies that `renv` does not automatically capture](../starting-your-analysis/determining-software-requirements/#pinning-dependencies-that-are-not-captured-automatically)
- These additional files and folders manage the `renv` environment, and you should not directly edit them:
    - `renv.lock`
    - The `renv` folder
    - `.Rprofile`


### Flags to create a Python module

#### `--use-jupyter`

!!! info
    The `hello-python` example module was created with this flag.
    Learn more about [using conda to manage your Python environment](../starting-your-analysis/determining-software-requirements/#managing-software-dependencies-in-python-with-conda).

- This flag will both:
    - Add a template Jupyter notebook to your module
      - To instead add a template Python script, use `--use-python`
    - Initialize a conda environment for your module
      - The conda environment will be named `openscpca-<module name>`
        - For example, if you name your module `celltype-ewings`, its conda environment will be named `openscpca-celltype-ewings`
      - The conda environment will include an installation of Jupyter that you will be able to launch with the `jupyter lab` command when the environment is active.

You can use this flag as:

```bash
./create-analysis-module.py my-module-name --use-jupyter

# Or, to create a Python script instead of a Jupyter notebook:
./create-analysis-module.py my-module-name --use-python
```

This will create skeleton module called `my-module-name` with these files and folders:

```
├── scripts
│   └── ...
├── results
│   └── README.md
├── plots
│   └── ...
├── scratch
│   └── ...
├── README.md
├── environment.yml
├── notebook-template.ipynb  # Or, script-template.py if you use `--use-python`
├── Dockerfile
├── .gitignore
└── .dockerignore
```

- You can use `notebook-template.ipynb` (or `script-template.py`) as a starting point for any Jupyter notebooks (or Python scripts) you create while writing your analysis
- You can use the `environment.yml` file to add packages to your module's conda environment



#### `--use-conda`

This flag will initialize a conda environment for your module, but will not add a template script or notebook.

You can use this flag as:

```bash
./create-analysis-module.py my-module-name --use-conda
```

This will create skeleton module called `my-module-name` with these files and folders:


```
├── scripts
│   └── ...
├── results
│   └── README.md
├── plots
│   └── ...
├── scratch
│   └── ...
├── README.md
├── environment.yml
├── Dockerfile
├── .gitignore
└── .dockerignore
```

- You can use the `environment.yml` file to add packages to your module's conda environment
