This file contains information about contributing to docs development.

`OpenScPCA-analysis` documentation is written in markdown with [`Material for MkDocs`](https://squidfunk.github.io/mkdocs-material/).
The [`mkdocs.yml`](../mkdocs.yml) is the documentation configuration file, and documentation pages are stored in this directory (`docs`).
To locally test documentation, you will need `python3` and dependencies listed in the [`requirements.txt`](./requirements.txt) file.
To install these dependencies, use your preferred approach:

```sh
# without conda
pip install -r requirements.txt

# with conda
conda install --file requirements.txt
```

To serve the documentation during local development, run `mkdocs serve` from the root directory of this repository.


## Docs organization

Documentation is written as a series of markdown files nested by topic in directories.

- All directories at the top-level inside `docs` will represent navbar sections.
  - Each top-level directory should contain an `index.md` file with an overall description of what that section contains.
  - The `index.md` should have an L1 header with the same title as the navbar section.
  - Markdown files in each directory represent sections shown along the left sidebar.
- Nested directories within each navbar section should be used to add an additional bold (gray) header on the left sidebar.
  - Markdown files within each nested directory will fall under this additional bold (gray) header.
- Any visual aids used in the docs should be placed in `docs/img`.

## Adding new documentation

_Instructions forthcoming on how to find the right place to add your docs._

Files and directories should be named with all lowercase letters using `-`, not `_`, where needed.
**Be sure to consult the [style guide](./general-style-guide.md) when writing your documentation.**
