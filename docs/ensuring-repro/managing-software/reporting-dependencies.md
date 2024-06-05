
## Reporting dependencies

All modules should have [clear documentation](../../contributing-to-analyses/analysis-modules/documenting-analysis.md) explaining how module dependencies are managed and how to activate any software environment(s).

In addition to this documentation, we recommend also recording your software dependencies directly from your code using some of the strategies detailed below.

## Reporting dependencies from R code

### Using `sessionInfo()` in R

The R programming language includes a helpful function called [`sessionInfo()`](https://stat.ethz.ch/R-manual/R-patched/library/utils/html/sessionInfo.html) that prints information about the R version and loaded packages.

#### In R Notebooks

Include `sessionInfo()` in the final chunk or cell of any computational notebook using R to ensure that the report output will include information about R and package versions.

#### In R scripts

Not using a notebook?
No problem!
You can capture the output of `sessionInfo()` in a text file by including the following in your script:

```r
sink("sessionInfo.txt")
sessionInfo()
sink()
```

This will result in the contents being written to a file called `sessionInfo.txt` in your working directory, so be sure to use unique filenames if multiple scripts in your analysis are writing out the contents of `sessionInfo()`.

You can post the contents of these files on [your pull requests](../../contributing-to-analyses/creating-pull-requests/index.md) or commit them to the repository.

## Reporting dependencies from Python code

### Using `session_info.show()` in Python

The [`session-info`](https://pypi.org/project/session-info/) Python package can be used to report version information about Python and loaded modules.
If you [created a module using `--use-jupyter` or `--use-python`](../../contributing-to-analyses/analysis-modules/creating-a-module.md#use-Jupyter), `session-info` was automatically included in the module's conda environment.

Import `session-info` by placing the following in [the Setup section of your Jupyter notebook](../../contributing-to-analyses/analysis-modules/notebook-structure.md#jupyter-notebooks) or in the [load packages section of your script](../../contributing-to-analyses/analysis-modules/script-structure.md):

```python
import session_info
```

#### In Jupyter notebooks

Include the following in the final cell of any Jupyter notebook using Python to ensure that the notebook reports version information:

```python
session_info.show()
```

#### In Python scripts

Using a script instead?

Import `contextlib` by placing the following in the [load packages section of your script](../../contributing-to-analyses/analysis-modules/script-structure.md):

```python
import contextlib
```

Then, you can use the following to write the output of `session_info.show()` to a text file called `session_info.txt` in your current directory:

```python
session_info_path = "session_info.txt"

with open(session_into_path, "w") as f:
    with contextlib.redirect_stdout(f):  # direct the session_info output to a file
        session_info.show(dependencies=True)
```

If you're using this approach with multiple scripts in your module, be sure to use unique file names when naming your session info text files.

You can post the contents of these files on [your pull requests](../../contributing-to-analyses/creating-pull-requests/index.md) or commit them to the repository.
