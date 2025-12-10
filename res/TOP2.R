library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(rrBLUP)
library(MASS)


# Read genotype 
Genotype<-read.table('00_Genotype.txt', header = TRUE, row.names = 1, sep = "\t")

# Read phenotype
Phenotype<-read.table('00_Phenotype.txt', header = TRUE, sep = "\t")

# Read environmental data
EnvFactors<-read.table('00_Envs.txt', header = TRUE, sep = "\t")

# Read known gene information
Gene_SNP_pos<-read.table('00_Gene_SNP_Pos.txt', header = T, sep = '\t')

Traits = unique(Phenotype$Traits)
TrainID = unique(Phenotype$line_code)
TestID = setdiff(colnames(Genotype), TrainID)

totalEnvs = c('A', 'B', 'C', 'D')
Envs = 'PreEnv'


############ Single environment prediction test set phenotype + construction of pseudo-true training set ############ 
source('SingleEnv_Pre.R')

Kinship<-A.mat(t(Genotype))
Kinship_Train<-Kinship[TrainID, TrainID]

Phe_Train_Pseudo<-data.frame()
Phe_Test<-data.frame()
for (trait in Traits) {
  res<-predict_BLUP_one_trait(trait, Genotype, Phenotype, Kinship, Kinship_Train)
  
  Phe_Train_Pseudo<-rbind(Phe_Train_Pseudo, res$Pseudo_Train)
  Phe_Test<-rbind(Phe_Test, res$Prediction_Test)
}


############ Predicting phenotypes across environments ############ 
source('NewEnv_Pre_Linear.R')
source('NewEnv_Pre_Nonlinear.R')

Phe_Train<-Phenotype %>% 
  group_by(Traits, Envs) %>% 
  mutate(Value = ifelse(is.na(Value), mean(Value, na.rm = TRUE), Value)) %>%
  ungroup() %>% 
  rename(Predicted = Value)

Phe_Envs_Train<-data.frame()
Phe_Envs_Train_Pseudo<-data.frame()
Phe_Envs_Test<-data.frame()
for (trait in Traits) {
  if (trait %in% c('DTA', 'DTS', 'DTT')) {
    # Training set
    phe<-predict_linear(trait, Phe_Train)
    Phe_Envs_Train<-rbind(Phe_Envs_Train, phe)
    # Pseudo-real training set
    phe<-predict_linear(trait, Phe_Train_Pseudo)
    Phe_Envs_Train_Pseudo<-rbind(Phe_Envs_Train_Pseudo, phe)
    # Testing set
    phe<-predict_linear(trait, Phe_Test)
    Phe_Envs_Test<-rbind(Phe_Envs_Test, phe)
  }else{
    # Training set
    phe<-predict_nonlinear(trait, Phe_Train)
    Phe_Envs_Train<-rbind(Phe_Envs_Train, phe)
    # Pseudo-real training set
    phe<-predict_nonlinear(trait, Phe_Train_Pseudo)
    Phe_Envs_Train_Pseudo<-rbind(Phe_Envs_Train_Pseudo, phe)
    # Testing set
    phe<-predict_nonlinear(trait, Phe_Test)
    Phe_Envs_Test<-rbind(Phe_Envs_Test, phe)
  }
}


############ Select materials ############ 
source('SmartTOP_Pre.R')
source('SmartGP.R')
source('ChooseMaterial.R')

trait = 'EW' # Target traits
trait_related = c('ERN', 'KNPR', 'KNPE', 'KWPE') # Traits that are strongly correlated with the target trait
selection_ratio = 0.02 # Identify ratio of individuals in the test set with maximized similarity to the target
improve_ratio = 0.1 # Change ratio of a trait value from the target
targetID = 'Test4718' # Any material can be specified as a target

MaterialID<-choose_material(trait, trait_related, selection_ratio, improve_ratio, targetID)
write.table(MaterialID, file = 'Choose_MaterialID.txt', quote = F, col.name = F, row.names = F, sep = "\t")
