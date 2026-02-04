# PR function
PR<-function(dt,A,nut){
  dt_M<-dt
  dt_M[nut]<-NA
  data<-data.frame(y = dt_M, gid = rownames(A))
  ans<-kin.blup(data = data, geno = "gid", pheno = "y", K = A)
  pl_te<-ans$g[nut] + mean(dt[-nut])
  Model = list(pcor_te = 0, pl_te = pl_te) 
  return(Model)
}


predict_BLUP_one_trait<-function(trait, Genotype, Phenotype, Kinship, Kinship_Train, Num_Fold = 10){
  set.seed(923)
  
  # Fill missing values
  Phenotype<-Phenotype %>% 
    group_by(Traits, Envs) %>% 
    mutate(Value = ifelse(is.na(Value), mean(Value, na.rm = TRUE), Value)) %>%
    ungroup()
  
  Pseudo_Train<-data.frame()
  Prediction_Test<-data.frame()
  Envs<-unique(Phenotype$Envs)
  
  for (env in Envs) {
    phe_sub<-Phenotype %>% filter(Traits == trait, Envs == env)
    train_vec<-phe_sub %>% filter(line_code %in% TrainID) %>% select(line_code, Value)
    
    GroupNum<-round(nrow(train_vec) / Num_Fold)

    ## Construct a pseudo-training set based on the specified number of folds.
    train_pred_all<-data.frame()
    for (i in 1:Num_Fold) {
      if (i < Num_Fold) {
        idx_test<-((i - 1) * GroupNum + 1):(i * GroupNum)
      } else {
        idx_test<-((i - 1) * GroupNum + 1):nrow(train_vec)
      }
      
      test_lines = train_vec$line_code[idx_test]
      train_lines = setdiff(train_vec$line_code, test_lines)
      
      Order = c(test_lines, train_lines)
      KA<-Kinship_Train[Order, Order]
      y<-c(rep(NA, length(test_lines)), train_vec$Value[match(train_lines, train_vec$line_code)])
      
      pred<-PR(dt = y, A = KA, nut = 1:length(test_lines))$pl_te
      pred_df <- data.frame(line_code = test_lines, Traits = trait, Envs = env, Predicted = pred)
      train_pred_all <- rbind(train_pred_all, pred_df)
    }
    
    Pseudo_Train<-rbind(Pseudo_Train, train_pred_all)
    
    ## Prediction testing set
    Order = c(TestID, TrainID)
    KA<-Kinship[Order, Order]
    train_value_df<-data.frame(line_code = as.character(TrainID)) %>% 
      left_join(train_vec, by = "line_code")
    y<-c(rep(NA, length(TestID)), train_value_df$Value)
    
    pred<-PR(dt = y, A = KA, nut = 1:length(TestID))$pl_te
    test_df<-data.frame(line_code = TestID, Traits = trait, Envs = env, Predicted = pred)
    Prediction_Test<-rbind(Prediction_Test, test_df)
  }
  
  return(list(Pseudo_Train = Pseudo_Train, Prediction_Test = Prediction_Test))
}

