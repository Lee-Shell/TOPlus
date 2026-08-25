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
**1. Genotype and environment independence** 
> Genotypes and environments are separated to prevent overlap in both genetic background and environmental conditions.  

**2. Model training for genotype performance**  
> Trait-specific models are trained using the training environments to estimate the performance of previously unobserved genotypes.

**3. Environment-driven factor selection and cross-environment prediction**  
> Training environments are iteratively treated as selection environments to identify key environmental factors associated with target traits. The remaining environments are used to model trait responses, allowing phenotype estimation in untested environments through aggregated environmental contributions.

**4. Target-oriented material prioritization**  
> Functional gene knowledge is integrated to estimate gene and trait weights, which are combined with predicted phenotypes to rank candidate hybrids relative to elite benchmarks under specific target environments.

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
