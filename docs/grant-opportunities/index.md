# Call for contributions to OpenScPCA and Terms for Grant Offering

As part of the OpenScPCA project, the Alex's Lemonade Stand Foundation (ALSF) Childhood Cancer Data Lab aims to characterize the ScPCA data, including assigning cell types to all samples across all tumor types.

We seek pediatric cancer researchers at non-profit institutions with experience analyzing single-cell RNA-seq datasets to help annotate and assign cell types to existing ScPCA datasets.

Contributors to OpenScPCA who annotate cell types:

-   May be eligible for a small one-time grant of research funds to be paid to their institution for advancing pediatric cancer research (see Award section below).

-   Will help create a set of publicly available cell type reference datasets across over 50 pediatric tumor types that the pediatric cancer research community can use to assign and label cell types in their own datasets.

-   Will have access to computational resources to complete the analyses.

For more information on the project, please refer to the [OpenScPCA documentation](https://openscpca.readthedocs.io/en/latest/).

All participation in this grant opportunity and contributions to the OpenScPCA project are subject to the Terms of Use, available at [https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/docs/policies/terms-of-use.md](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/docs/policies/terms-of-use.md), and all Policies posted at [https://openscpca.readthedocs.io/en/latest/policies/](https://openscpca.readthedocs.io/en/latest/policies/). Each researcher and the institution to which the researcher is affiliated must agree to all of the terms of the Terms of Use and all Policies, which are incorporated herein by reference as though fully set forth herein.

## Cell type annotation

We seek pediatric cancer researchers at non-profit institutions to help annotate cell types for all samples in a given group as identified by the sample table. We expect final cell type annotations will include tumor (malignant) and normal (non-malignant) classifications and specific cell type labels for non-malignant cells (e.g., T cell, fibroblast) and cell states present in tumor cells (e.g., mesenchymal-like cells).

## Award

### Available Grants

A maximum of 27 grants are available. A complete table of available grants is available at [https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/571](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/571).

This table contains a unique group identifier, project identifier, project title, list of diagnoses, and total number of samples for each group. Groups of samples encompass samples from a single ScPCA project representing either a single diagnosis or a group of related diagnoses.

The table includes the maximum amount that may be awarded in USD for annotating all samples in each group, with the award amount based on the total number of samples in that group as follows:

-   $4000 for annotating groups with 30 or more samples

-   $2000 for annotating groups with 16 - 29 samples

-   $1000 for annotating groups with 15 or fewer samples

Only one grant will be awarded per group of samples and will be awarded to the institution for the first researcher to complete cell type annotation satisfactorily and meets all eligibility and other requirements established by ALSF.

ALSF will update the table of available grants at [https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/571](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/571) to indicate which grant awards have been awarded and are no longer available.

ALSF will only award grants for submissions meeting the eligibility requirements and acceptance criteria. All decisions regarding applications and grants, including satisfaction of requirements, are made solely by ALSF, in its sole discretion.

### Distribution and Permitted Use of Funds

-   The funds will be distributed to the researcher's institution and must be used for pediatric cancer research. If the institution is not conducting such research at the time the grant is to be distributed, the grant will be deemed declined.

-   ALSF does not allow the use of funds for research utilizing human embryonic stem cells or nonhuman primates. Research with human induced pluripotent stem cells is permissible.

-   The institution receiving the grant is responsible for compliance with any governmental reporting requirements and related tax obligations, if any.

### General Eligibility Requirements

-   Researchers must be employees of or have an appointment at a non-profit institution when they fill out the award application. Grants will be paid to the institution identified in the grant award application.

-   Researchers with a primary affiliation at an institution in a geographic location where United States-based institutions cannot send funds due to governmental policy are ineligible. For more information, please see the current U.S. Department of the Treasury Sanction Programs and Country information at [https://ofac.treasury.gov/sanctions-programs-and-country-information](https://ofac.treasury.gov/sanctions-programs-and-country-information).

-   An individual researcher can qualify for a maximum of three grants.

-   There is no limit to the number of grants that will be awarded to an individual institution.

### Submission Acceptance Criteria

To be eligible for a grant to the researcher's institution, researchers must complete cell type annotation for all samples assigned to a specified group identifier.

Researchers must complete the cell type annotation analysis openly and collaboratively by contributing their analysis to the OpenScPCA project. The analysis code and table described below must be available in the public [OpenScPCA-analysis repository](https://github.com/AlexsLemonade/OpenScPCA-analysis).

Annotations must be scientifically justified (i.e., data-driven or derived using marker genes supported by scientific literature). All annotations, underlying code, and results will be evaluated for scientific justification, correctness, and integrity through a peer review process conducted openly on GitHub. An ALSF Childhood Cancer Data Lab staff member must approve the researchers' contributions, and all decisions are in the sole discretion of ALSF.

Researchers must submit a table in comma-separated value (CSV) or tab-separated value (TSV) format containing the final cell type labels for all samples in a group with the following columns:

-   `scpca_sample_id`: The unique identifier for the sample, as provided on the ScPCA Portal.

-   `cell_barcode`: 10x Genomics cell barcode associated with each cell.

-   `tumor_cell_classification`: Classification of each cell as either tumor (malignant) or normal (non-malignant). The values of this column should be either "tumor" or "normal."

-   `cell_type_assignment`: The cell type or label assigned to the cell, including the non-malignant cell types (T cell, fibroblast, etc.) and any cell states present in tumor cells (e.g., mesenchymal-like cells).

-   `CL_ontology_id`: (*Optional*) The [cell ontology (CL) identifier](https://www.ebi.ac.uk/ols4/ontologies/cl) associated with the designated cell type. Although these are not required, they are strongly encouraged to help standardize cell type assignment across datasets available on the Portal.

If researchers use marker genes during the process of annotating cell types, a table with the following columns is also required for final submission:

-   `ensembl_gene_id`: Ensembl gene identifier for the marker gene.

-   `cell_type`: Cell type label the marker gene is associated with. This value must match the value in the `cell_type_assignment` column in the above annotations table.

-   `Source`: The source from which the marker genes were obtained. Typically, this will be a DOI for a publication using this marker gene.

### Important Dates

Researchers may submit tables of annotations starting at 9 AM Eastern Daylight Time on July 1, 2024.

The final table must be approved via pull request by **5 PM Eastern Daylight Time on November 1, 2024.**

Researchers should plan for multiple rounds of peer review before the final table submission. Please start participating as early as possible to avoid missing the submission deadline.

### Notification of Eligibility and Award

Researchers who have completed cell type annotation for a group of samples and are eligible for a grant will be notified via email using the email address supplied on [the OpenScPCA interest form](https://share.hsforms.com/1MlLtkGYSQa6j23HY_0fKaw336z0).

**Researchers must complete the award application form to notify ALSF of their acceptance of the award within ten days of notification of eligibility or by 5 PM Eastern Standard Time on November 15, 2024, whichever is sooner.**

Upon application approval, researchers will receive an award notification letter via email at the email address supplied on the OpenScPCA interest form.

### Contact

Questions about the project, grant requirements, or grant awards may be directed to the ALSF Childhood Cancer Data Lab team at [scpca@ccdatalab.org](mailto:scpca@ccdatalab.org).

##  

## Getting Started

If you want to annotate one of the groups of samples indicated by the group identifier on the samples table, please let us know by filing a [GitHub discussion to propose a new analysis](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/new?category=propose-a-new-analysis).

In your post, please mention the unique group identifier for which you would like to add cell type annotations along with the project title, project identifier, and list of diagnoses for that group. You should also include an outline of the approach you plan to take to complete the cell type annotation, including an explanation of the methods you plan to use to perform cell type annotation. An ALSF Childhood Cancer Data Lab staff member will then respond to your post and help you set up to start your analysis and contribution to the OpenScPCA project.
