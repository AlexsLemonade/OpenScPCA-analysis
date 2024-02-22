This file contains information about contributing docs development.

`OpenScPCA-analysis` documentation is written in markdown with [`Material for MkDocs`](https://squidfunk.github.io/mkdocs-material/).
The [`mkdocs.yml`](../mkdocs.yml) is the documentation configuration file, and documentation pages are stored in the `docs` directory.
To locally test documentation, you will need `python3` and dependencies listed in the [`requirements.txt`](./requirements.txt) file.
To install these dependencies, use your preferred approach:

```sh
# without conda
pip install -r requirements.txt

# with conda
conda install --file requirements.txt
```


To serve the documentation, run `mkdocs serve` from the root directory of this repository.
