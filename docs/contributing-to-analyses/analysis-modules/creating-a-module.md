# Creating an analysis module

Once you have [discussed your proposed analysis](../../communications-tools/index.md#github-discussions) and [filed an issue to get started](../../communications-tools/index.md#github-issues), you're ready to create an analysis module.

We have provided a script [`create-analysis-module.py`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/create-analysis-module.py) you can use to make this process easier.
Before you begin, please review [this documentation](../analysis-modules/index.md#skeleton-analysis-module-contents) for an overview of the skeleton file structure of an empty analysis module.


!!! note "Is this your first analysis module?"
    If you are creating your first analysis module, we recommend that you [scope your first pull request](../creating-pull-requests/scoping-pull-requests.md) to contain only:

    - The analysis module skeleton created by running `create-analysis-module.py`
    - Some [initial documentation of your module](./documenting-analysis.md) in the `README.md` file

    This way, you can get some experience filing pull requests, undergoing code review, and using other Git skills before you really get going on the analysis.


## Running the module creation script

You should use the [`create-analysis-module.py`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/create-analysis-module.py) script to create a skeleton module for your analysis in the `analyses` folder.

Before running this script, you should determine two things:

- The name of your module
    - You will provide this name as an argument to the script
-  Whether you will use R or Python (or both!) in your module
    - You can use certain script flags to add a template notebook and/or set up an environment (`renv` or conda) for your chosen language, [as described below in detail](#module-creation-script-flags)
        - If you do not use additional flags, the script will create a [skeleton module](../analysis-modules/index.md#skeleton-analysis-module-contents) only


To run this script, take these steps:

1. Open a [terminal](../../software-platforms/general-tools/using-the-terminal.md) window.
    - You may wish to launch terminal [from GitKraken](../../software-platforms/general-tools/using-the-terminal.md#gitkraken) so that you are automatically placed in the repository.

2. Make sure you are working in your `base` conda environment by running `conda activate base`.

3. As needed, navigate using `cd` to your `OpenScPCA-analysis` fork.

4. Run the module creation script:

    ```bash
    # Create a module called `my-module-name`
    ./create-analysis-module.py my-module-name <additional flags here>
    ```

      - The script requires one first argument: the name of your module.
          - The example code above will creates an analysis module called `my-module-name`
      - You can also use an [additional flag](#module-creation-script-flags) to add a template notebook or set up an environment with `renv` or conda


!!! tip
    You can always run the script with the `--help` flag to reveal the help menu and see other flags you can use:

    ```bash
    ./create-analysis-module.py --help
    ```

## Module creation script flags


We recommend using one of these flags when creating your module.
Each flag will create a skeleton module with the given files and folders.

Note that you can use both an R and a Python flag if you want to write your module in both languages.

Each example below shows the resulting module file structure when using each flag.

### Flags to create an R module

#### The `--use-r` flag

Use this flag to add a template R notebook to your module:

```bash
# Create a module called `my-module-name` with a template R Markdown notebook
./create-analysis-module.py my-module-name --use-r
```

<div class="grid" markdown>

- You can use `notebook-template.Rmd` as a starting point for any R Markdown notebooks you create while writing your analysis


```{ .console .no-copy title="Module directory with --use-r flag"}
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




#### The `--use-renv` flag

!!! info
    The `hello-r` example module was created with this flag.
    Learn more about [using `renv` to manage your R environment](../determining-requirements/determining-software-requirements.md#using-renv).

Use this flag to:

- Add a template R notebook to your module
- Initialize an `renv` environment for your module


```bash
# Create a module called `my-module-name` with a template R Markdown notebook
#   and an `renv` environment
./create-analysis-module.py my-module-name --use-renv
```

<div class="grid" markdown>

- You can use `notebook-template.Rmd` as a starting point for any R Markdown notebooks you create while writing your analysis
- You can use `components/dependencies.R` to [pin R package dependencies that `renv` does not automatically capture](../determining-requirements/determining-software-requirements.md#pinning-dependencies-that-are-not-captured-automatically)
- These additional files and folders manage the `renv` environment, and you should not directly edit them:
    - `renv.lock`
    - The `renv` folder
    - `.Rprofile`


```{ .console .no-copy title="Module directory with --use-renv flag"}
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

</div>

### Flags to create a Python module

#### The `--use-jupyter` flag

!!! info
    The `hello-python` example module was created with this flag.
    Learn more about [using conda to manage your Python environment](../determining-requirements/determining-software-requirements.md#managing-software-dependencies-in-python-with-conda).

Use this flag to:

- Add a template Jupyter notebook to your module
    - To add a template Python script instead, use `--use-python`
- Initialize a conda environment for your module
    - The conda environment will be named `openscpca-<module name>`
    - For example, if you name your module `celltype-ewings`, its conda environment will be named `openscpca-celltype-ewings`
        - The conda environment will include an installation of Jupyter that you can launch with the `jupyter lab` command when the environment is active

```bash
# Create a module called `my-module-name` with a template Jupyter notebook
#   and a conda environment with Jupyter installed
./create-analysis-module.py my-module-name --use-jupyter

# Or, create a module called `my-module-name` with a template python script
#   and a conda environment (Jupyter not installed)
./create-analysis-module.py my-module-name --use-python
```

<div class="grid" markdown>


- You can use `notebook-template.ipynb` (or `script-template.py`) as a starting point for any Jupyter notebooks (or Python scripts) you create while writing your analysis
- You can use the `environment.yml` file to add additional packages to your module's conda environment


```{ .console .no-copy title="Module directory with --use-jupyter flag"}
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

</div>

#### The `--use-conda` flag

Use this flag to initialize a conda environment in your module, but without a template script or notebook.
The conda environment will be named `openscpca-<module name>`.


```bash
# Or, create a module called `my-module-name` with a conda environment
./create-analysis-module.py my-module-name --use-conda
```

<div class="grid" markdown>

- You can use the `environment.yml` file to add packages to your module's conda environment


```{ .console .no-copy title="Module directory with --use-conda flag"}
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

</div>