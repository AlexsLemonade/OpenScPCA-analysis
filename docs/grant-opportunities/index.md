# Call for Contributions to OpenScPCA and Terms for Grant Offering

As part of the OpenScPCA project, the Alex's Lemonade Stand Foundation (ALSF) Childhood Cancer Data Lab aims to characterize the ScPCA data, including assigning cell types to all samples across all tumor types.

We seek pediatric cancer researchers at non-profit institutions with experience analyzing single-cell RNA-seq datasets to help annotate and assign cell types to existing ScPCA datasets.

Contributors to OpenScPCA who annotate cell types:

-   May be eligible for a small one-time grant of research funds to be paid to their institution for advancing pediatric cancer research (see Award section below).

-   Will help create a set of publicly available cell type reference datasets across over 50 pediatric tumor types that the pediatric cancer research community can use to assign and label cell types in their own datasets.

-   Will have access to computational resources to complete the analyses.

For more information on the project, please refer to the [OpenScPCA documentation](https://openscpca.readthedocs.io/en/latest/).

All participation in this grant opportunity and contributions to the OpenScPCA project are subject to the Terms of Use, available at <https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/docs/policies/terms-of-use.md>, and all Policies posted at <https://openscpca.readthedocs.io/en/latest/policies/>. 
Each researcher and the institution to which the researcher is affiliated must agree to all of the terms of the Terms of Use and all Policies, which are incorporated herein by reference as though fully set forth herein.

## Cell Type Annotation

We seek pediatric cancer researchers at non-profit institutions to help annotate cell types for all samples in a given group as identified by the sample table.
We expect final cell type annotations will include tumor (malignant) and normal (non-malignant) classifications and specific cell type labels for non-malignant cells (e.g., T cell, fibroblast) and cell states present in tumor cells (e.g., mesenchymal-like cells).

## Awards

### Available Grants

A maximum of 19 grants are available. 
Each grant is tied to a specific sample group assigned a unique group identifier. 
A complete table of group identifiers is available at <https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/1062>.

This table contains the unique group identifier, project identifier, project title, list of diagnoses, and total number of samples for each group. 
Groups of samples encompass samples from a single ScPCA project representing either a single diagnosis or a group of related diagnoses.

Grant opportunities are available in two categories with distinct submission acceptance requirements: the Advisor Award and the Analyst Award.

For an individual group, the award categories are mutually exclusive.
Only one grant will be awarded per group of samples. 
It will be awarded to the institution for the first researcher to complete the submission acceptance requirements satisfactorily for either the Advisor or Analyst Award and meet all eligibility and other requirements established by ALSF.

ALSF will update the table of available grants at <https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/1062> to indicate which grant awards have been awarded and are no longer available.

ALSF will only award grants for submissions meeting the eligibility and acceptance criteria. 
All decisions regarding applications and grants, including satisfaction of requirements, are made solely by ALSF, in its sole discretion.

#### The Analyst Award 

The maximum amount that may be awarded in USD for an Analyst Award is determined by the number of samples in a group as follows:

-   $4000 for annotating groups with 30 or more samples

-   $2000 for annotating groups with 16 - 29 samples

-   $1000 for annotating groups with 15 or fewer samples

The table of available grants at <https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/1062> indicates the maximum amount that may be awarded for an Analyst Award for a given group.

#### The Advisor Award

The maximum amount that may be awarded for an Advisor Award is $500 USD. 
The table of available grants at <https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/1062> indicates the maximum amount that may be awarded for an Advisor Award for a given group.

### Distribution and Permitted Use of Funds

-   The funds will be distributed to the researcher's institution and must be used for pediatric cancer research. 
    If the institution is not conducting such research at the time the grant is to be distributed, the grant will be deemed declined.

-   ALSF does not allow the use of funds for research utilizing human embryonic stem cells or nonhuman primates. 
    Research with human induced pluripotent stem cells is permissible.

-   The institution receiving the grant is responsible for compliance with governmental reporting requirements and related tax obligations, if any.

### General Eligibility Requirements

-   Researchers must be employees of or have an appointment at a non-profit institution when they fill out the award application.
    Grants will be paid to the institution identified in the grant award application.

-   Researchers with a primary affiliation at an institution in a geographic location where United States-based institutions cannot send funds due to governmental policy are ineligible. 
    For more information, please see the current U.S. Department of the Treasury Sanction Programs and Country information at <https://ofac.treasury.gov/sanctions-programs-and-country-information>.

-   An individual researcher can qualify for a maximum of three grants.

-   There is no limit to the number of grants that will be awarded to an individual institution.

### Submission Acceptance Criteria

#### The Analyst Award

To be eligible for an Analyst Award grant to the researcher's institution, researchers must complete cell type annotation for all samples assigned to a specified group identifier.

Researchers must complete the cell type annotation analysis openly and collaboratively by contributing their analysis to the OpenScPCA project.
The analysis code and table described below must be available in the public [OpenScPCA-analysis repository](https://github.com/AlexsLemonade/OpenScPCA-analysis).

Annotations must be scientifically justified (i.e., data-driven or derived using marker genes supported by scientific literature).
Samples in a group may be excluded from cell type annotation with justification (e.g., the sample's library consists of very few viable cells) and only with prior approval from ALSF.

All annotations, underlying code, and results will be evaluated for scientific justification, correctness, and integrity through a peer review process conducted openly on GitHub.
An ALSF Childhood Cancer Data Lab staff member must approve the researchers' contributions, and all decisions are in the sole discretion of ALSF.

Please be aware that multiple automated cell type annotation methods have already been applied to ScPCA data. 
If an automated cell type annotation method is used in the analysis, researchers will be asked to demonstrate its superiority over the existing labels during peer review.
For example, the following characteristics may demonstrate superior annotation: greater specificity in labeled cell types, fewer cells labeled as unknown, and/or more coherent expression of curated marker genes within labeled cell types.

Researchers should plan for multiple rounds of peer review before their submission is accepted. 
Please start participating as early as possible to avoid missing the submission deadline in [Important Dates](#important-dates).

Researchers must submit a table in comma-separated value (CSV) or tab-separated value (TSV) format containing the final cell type labels for all samples in a group with the following columns:

-   `scpca_sample_id`: The unique identifier for the sample, as provided on the ScPCA Portal.

-   `cell_barcode`: 10x Genomics cell barcode associated with each cell.

-   `tumor_cell_classification`: Classification of each cell as either tumor (malignant) or normal (non-malignant). 
    The values of this column should be either "tumor" or "normal."

-   `cell_type_assignment`: The cell type or label assigned to the cell, including the non-malignant cell types (T cell, fibroblast, etc.) and any cell states present in tumor cells (e.g., mesenchymal-like cells). 
    No more than 10% of the cells in a specified group may be labeled as unknown without explicit justification approved by ALSF.

-   `CL_ontology_id`: (*Optional*) The [cell ontology (CL) identifier](https://www.ebi.ac.uk/ols4/ontologies/cl) associated with the designated cell type.
    Although these are not required, they are strongly encouraged to help standardize cell type assignment across datasets available on the Portal.

If researchers use marker genes during the process of annotating cell types, a table with the following columns is also required for final submission:

-   `ensembl_gene_id`: Ensembl gene identifier for the marker gene.

-   `cell_type`: Cell type label the marker gene is associated with. 
    This value must match the value in the cell_type_assignment column in the above annotations table.

-   `Source`: The source from which the marker genes were obtained.
    Typically, this will be a DOI for a publication using this marker gene.

Researchers must also complete the Description section in the README of the analysis directory that contains the analytical code used to perform the cell type annotation.
The Description section must describe the steps taken to complete the cell type annotation, including the types of analyses and methods or tools used.

#### The Advisor Award

To be eligible for an Advisor Award grant to the researcher's institution, researchers must submit and answer queries regarding a detailed protocol for cell type annotation for all diagnoses assigned to a specified group identifier. 
In contrast to the Analyst Award grant, researchers qualifying for an Advisor Award are not expected to perform analyses themselves.

Researchers must post the protocol as a GitHub Discussions post in the [OpenScPCA-analysis repository](https://github.com/AlexsLemonade/OpenScPCA-analysis) using [the Submit a cell type annotation protocol template](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/new?category=submit-a-detailed-protocol).

Protocols must be scientifically justified and detailed enough for an ALSF Childhood Cancer Data Lab staff member to complete cell type annotation as described. 
Protocols that do not include citations of existing literature or publicly available analyses will not be accepted, and citations of personal communications are insufficient to qualify.

Researchers must satisfactorily answer any questions posted as comments on their GitHub Discussions post before their submission will be approved. 
Researchers should anticipate that there may be multiple rounds of questions about their protocol.
Please start participating as early as possible to avoid missing the submission deadline in [Important Dates](#important-dates).

As indicated by a GitHub Discussions comment, an ALSF Childhood Cancer Data Lab staff member must approve the researchers' contributions. 
All decisions are in the sole discretion of ALSF.

### Important Dates

Researchers may submit materials starting at 9 AM Eastern Standard Time on March 3, 2025.

*Analyst Awards*: The final table and Description section of the README must be approved via pull request by **5 PM Eastern Daylight Time on October 31, 2025.**

*Advisor Awards*: The final protocol must be approved, as indicated by a GitHub Discussions comment from an ALSF Childhood Cancer Data Lab staff member, by **5 PM Eastern Daylight Time on October 31, 2025.**

### Notification of Eligibility and Award

Researchers who are eligible for a grant will be notified via email using the email address supplied on [the OpenScPCA interest form](https://share.hsforms.com/1MlLtkGYSQa6j23HY_0fKaw336z0).

**Researchers must complete the award application form to notify ALSF of their acceptance of the award within ten days of notification of eligibility or by 5 PM Eastern Standard Time on November 14, 2025, whichever is sooner.**

Upon application approval, researchers will receive an award notification letter via email at the email address supplied on the OpenScPCA interest form.

### Contact

Questions about the project, grant requirements, or grant awards may be directed to the ALSF Childhood Cancer Data Lab team at <scpca@ccdatalab.org>.


## Getting Started 

#### The Analyst Award 

If you wish to annotate one of the groups of samples indicated by the group identifier on the samples table as part of an Analyst Award submission, please let us know by filing a [GitHub Discussion to propose a new analysis](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/new?category=propose-a-new-analysis).

In your post, please mention the unique group identifier for which you would like to add cell type annotations, along with the project title, project identifier, and list of diagnoses for that group.
You should also include an outline of the approach you plan to take to complete the cell type annotation, including an explanation of the methods you plan to use to perform cell type annotation.
An ALSF Childhood Cancer Data Lab staff member will then respond to your post and help you set up to start your analysis and contributions to the OpenScPCA project.

#### The Advisor Award

If you wish to submit a detailed protocol as part of an Advisor Award submission, please use [the Submit a cell type annotation protocol template](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/new?category=submit-a-detailed-protocol) to post on GitHub Discussions.

In your post, please mention the unique group identifier for which you are adding a protocol for cell type annotations, along with the project title, project identifier, and list of diagnoses for that group.

An ALSF Childhood Cancer Data Lab staff member will respond to your post with any questions or clarifications.
