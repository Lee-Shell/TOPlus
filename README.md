# TOPlus

Prioritizing the right hybrids for the right environments before large-scale field testing.

<img width="3515" height="1757" alt="TOPlus" src="https://github.com/user-attachments/assets/9d517e34-e68c-4ef6-9386-48e934374550" />

TOPlus is an environment-targeted hybrid prioritization framework that integrates enviromic data with functional gene knowledge to support multi-trait selection in defined target environments. By leveraging cross-environment signals and prior biological information, TOPlus enables reliable prediction of hybrid performance across untested genotypes and environments while identifying candidate materials that deliver environment-specific advantages without compromising overall trait balance. The framework provides a knowledge-guided and environment-aware foundation for predictive breeding and deployment decisions under complex genotype-by-environment scenarios.

## Key Features
**1. Environment-targeted prioritization**  
> Identifies candidate hybrids optimized for specific agroecological conditions rather than relying on broad, environment-agnostic selection.

**2. Knowledge-guided selection**    
> Incorporates functional gene priors to translate biological insight into actionable breeding decisions.

**3. Cross-environment inference**  
> Leverages environmental correlations to support prediction and deployment in previously untested environments.

**4. Balanced multi-trait optimization**  
> Selects materials that achieve environment-specific gains while maintaining overall agronomic stability.

**5. Deployment-oriented framework**   
> Supports region-level suitability assessment and extrapolative placement of elite hybrids.

## Why TOPlus?
Modern breeding generates abundant genomic discoveries, yet translating these findings into deployment-ready cultivar decisions remains a major bottleneck. TOPlus addresses this gap by unifying genomic signals, environmental heterogeneity, and biological knowledge within a single prioritization framework — enabling breeders to move from prediction to actionable selection with greater confidence.

## Methodological Workflow
TOPlus operates through four major stages:  

**1. Phenotype prediction**  
> Genomic and environment-specific phenotypic information from the training population is first used to predict candidate performance. For prediction in untested environments, training environments are iteratively treated as selection environments to identify informative environmental variables, and predictions across selection rounds are aggregated to estimate phenotypic performance in the target environment.  

**2. Genotype fitting**  
> Prior functional-gene variants are evaluated in the training population, and fitted gene values are estimated using a mixed linear model based on multi-trait phenotypic information.  

**3. Weight learning**  
> Observed and predicted phenotypes, together with observed and fitted gene information, are used to learn environment-specific weights for multiple traits and functional-gene priors.  

**4. Target-oriented material prioritization**  
> The learned weights are used to calculate the integrated phenotypic-genotypic similarity between each candidate hybrid and the target hybrid. Candidates with high similarity to the target while showing superior predicted performance in the target trait are prioritized for breeding selection.

## Installation
TOPlus requires R (≥ 4.2.0) and the following packages:
```
install.packages(c("dplyr", "tidyr", "stringr", "tibble", "rrBLUP", "MASS"))
```

## Quick Start
Run the TOPlus pipeline via the master script:
```
source("TOPlus.R")
```
This script automatically executes all modules in the correct order.  
Executing this script reproduces the complete TOPlus prioritization workflow.

## Input Data

TOPlus requires four types of input data.

### 1. Genotype data
The genotype file is a marker-by-material matrix. Rows represent SNP markers and columns represent individual materials. Training and candidate materials should be included in the same matrix.

### 2. Phenotype data
Phenotypic records are provided in long format with four required columns:

- `line_code`: material identifier
- `Traits`: trait identifier
- `Envs`: environment identifier
- `Value`: observed phenotypic value

Materials present in the phenotype file are treated as the training population, whereas genotype columns without phenotype records are treated as candidate materials.

### 3. Environmental data
Environmental variables should be provided for all training environments and the target environment. Environment names must correspond to those specified in the analysis settings.

### 4. Functional-gene information
The functional-gene file provides genomic coordinates used to assign SNP markers to prior genes. Marker coordinates and gene coordinates must use the same reference genome.

## User-defined Settings

Before running TOPlus, users should specify the breeding scenario in the master script.

- `totalEnvs`: training environments
- `Envs`: target environment
- `trait`: primary target trait
- `trait_related`: additional traits considered during prioritization
- `targetID`: target material
- `selection_ratio`: proportion of candidate materials retained
- `improve_ratio`: desired proportional change in the target trait

These settings correspond to the demonstration workflow and should be modified according to the user's breeding objective and dataset.

## Pipeline Structure
TOPlus is implemented as a modular pipeline orchestrated by a master script.  
The pipeline automatically selects linear or nonlinear prediction strategies depending on trait architecture.  
```text
TOPlus/
├── TOPlus.R
├── Data_demo/
│   ├── 00_Genotype_Demo.txt
│   ├── 00_Phenotype_Demo.txt
│   ├── 00_Envirotyping_Demo.txt
│   └── 00_Gene_SNP_Pos_Demo.txt
├── res/
│   ├── SingleEnvironmentPrediction.R
│   ├── LinearPrediction.R
│   ├── NonlinearPrediction.R
│   └── MaterialPrioritization.R
├── README.md
└── LICENSE
```

## Output

After successful execution, TOPlus generates the final candidate-prioritization file:

`TOPlus_Recommended_Materials.txt`

This file contains the candidate materials prioritized by TOPlus under the user-defined target environment and breeding objective.

Intermediate prediction and model-fitting objects are generated during execution and can be inspected in the R session for more detailed analyses.

## Using Your Own Data

To apply TOPlus to a new dataset:

1. Replace the demonstration genotype, phenotype, environmental, and functional-gene files with user-defined datasets following the same format.
2. Ensure that material identifiers are consistent between genotype and phenotype data.
3. Ensure that environment identifiers are consistent between phenotype and environmental data.
4. Ensure that SNP positions and prior-gene coordinates use the same reference genome.
5. Modify the training environments, target environment, target trait, related traits, and target material in the master script.
6. Run the master script and inspect the generated candidate-ranking output.

## Data Availability

The demonstration datasets required to run the example workflow are provided in the `Data_demo` directory.

The genotype and phenotype datasets used in the study are available from Zenodo:
https://zenodo.org/records/18605163

## Citation

If you use TOPlus in your research, please cite:

Fang et al. (2026). TOPlus: An algorithm integrating prior genetic and environmental information for multi-trait synergic selection in maize hybrid breeding. *Plant Communications*. https://doi.org/10.1016/j.xplc.2026.102089
