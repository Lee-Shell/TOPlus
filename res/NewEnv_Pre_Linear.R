predict_linear<-function(trait, df_phe_train){
  # Phenotype Range Threshold
  P_range = 0.1
  # Fixed Threshold
  C_cutoff = 0.25
  
  Phe0<-df_phe_train %>%
    filter(Traits == trait) %>%
    dplyr::select(-Traits) %>%
    pivot_wider(names_from = Envs, values_from = Predicted) %>%
    arrange(as.numeric(str_extract(line_code, "\\d+"))) %>%
    column_to_rownames(var = "line_code")
  
  # Phenotype range Filtered
  range_max = ceiling(max(Phe0) * (1 + P_range))
  range_min = floor(min(Phe0) * (1 - P_range))
  
  factors0<-EnvFactors %>%
    filter(Traits == trait) %>%
    column_to_rownames(var = colnames(.)[1]) %>%
    dplyr::select(all_of(c(totalEnvs, Envs)))

  df_phe_PreEnv<-data.frame()
  for (testEnvs in totalEnvs) {
    df_predict_testEnvs<-matrix(NA, nrow = nrow(Phe0), ncol = nrow(factors0),
                                dimnames = list(rownames(Phe0), rownames(factors0)))
    df_predict_Envs<-df_predict_testEnvs
    
    for (FactorID in rownames(factors0)) {
      df_Factor<-t(factors0[FactorID, , drop = FALSE])
      colnames(df_Factor) = 'Factor'
      for (MaterID in rownames(Phe0)) {
        # merge
        df_Mater<-t(Phe0[MaterID, , drop = FALSE])
        colnames(df_Mater) = 'Mater'
        
        df_mid<-merge(df_Factor, df_Mater, by = "row.names", all = TRUE)
        rownames(df_mid) = df_mid[, 1]
        df_mid<-df_mid[, -1]
        
        # Fitting environments other than testEnvs and Envs
        df_test<-df_mid[-which(rownames(df_mid) %in% c(testEnvs, Envs)), ]
        model<-lm(Mater ~ Factor, data = df_test)
        # Predicting phenotypes of testEnvs and Envs
        df_predict_testEnvs[MaterID, FactorID] = predict(model, 
                                                         newdata = data.frame(x = df_mid[testEnvs, 'Factor', drop = FALSE]))
        df_predict_Envs[MaterID, FactorID] = predict(model, 
                                                     newdata = data.frame(x = df_mid[Envs, 'Factor', drop = FALSE]))
      }
    }
    
    # Keep factors whose population prediction values are within the range
    df_predict_testEnvs<-as.data.frame(df_predict_testEnvs)
    df_predict_testEnvs<-df_predict_testEnvs[, sapply(df_predict_testEnvs,
                                                      function(x) all(x >= range_min & x <= range_max)), 
                                             drop = FALSE]
    if (ncol(df_predict_testEnvs) == 0) {
      next
    }
    
    # Add testEnvs true phenotype
    df_predict_testEnvs$Raw<-Phe0[, which(colnames(Phe0) %in% testEnvs)]
    # Cor
    cors<-sapply(df_predict_testEnvs[, -ncol(df_predict_testEnvs)], 
                 function(col) cor(col, df_predict_testEnvs$Raw))
    df_cor_testEnvs<-data.frame(FactorsRange = rownames(factors0)[seq_along(cors)], 
                                Cor = cors, row.names = NULL)
    # Keep factors that pass the threshold
    df_cor_testEnvs<-df_cor_testEnvs[df_cor_testEnvs$Cor >= C_cutoff, ,drop = FALSE]
    
    if (nrow(df_cor_testEnvs) > 0) {
      df_predict_Envs<-t(df_predict_Envs[, df_cor_testEnvs$FactorsRange, drop = FALSE])
      df_phe_PreEnv<-rbind(df_phe_PreEnv, df_predict_Envs)
    }
  }
  
  Phe_PreEnv<-data.frame(line_code = colnames(df_phe_PreEnv),
                         Predicted = colMeans(df_phe_PreEnv), 
                         Traits = trait)
  return(Phe_PreEnv)
}
