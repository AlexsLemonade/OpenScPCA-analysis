# Normal references

This directory stores RDS files with SCE objects to use as normal cell references with `inferCNV`, organized by project.
The RDS files are not stored in version control due to their large size.
Each project directory will also have a TSV tile listing all of the consensus cell types included in each reference, which is included in version control.

All files are created by project-specific scripts found in `../../scripts/build-normal-reference/`.

The TSV file `reference-cell-groups.tsv` contains consensus cell types organized into reference cell group categories, and it is used to help create references across projects.
The exception to this is `SCPCP000015`; this directory contains its own file `reference-celltypes.tsv` giving which cell types belong to each normal reference for this project.
