predict_nonlinear<-function(trait, df_phe_train){
  # Phenotype Range Threshold
  P_range = 0.1
  # Fixed Threshold
  C_cutoff = 0.25
  
  Phe0<-df_phe_train %>%
    filter(Traits == trait) %>%
    dplyr::select(-Traits) %>%
    pivot_wider(names_from = Envs, values_from = Predicted) %>%
    arrange(as.numeric(str_extract(line_code, "\\d+")))
  # Phenotype range Filtered
  range_max = ceiling(max(Phe0[, 2:5]) * (1 + P_range))
  range_min = floor(min(Phe0[, 2:5]) * (1 - P_range))
  
  # Keep the best R for each factor
  factors0<-EnvFactors %>%
    filter(Traits == trait) %>%
    group_by(FactorsRange) %>%
    slice_max(order_by = abs(R), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    column_to_rownames(var = colnames(.)[1])
  
  df_phe_PreEnv<-data.frame()
  for (testEnvs in totalEnvs) {
    factors_M<-factors0[, c(setdiff(totalEnvs, testEnvs), testEnvs)]
    Phe_M<-Phe0[, colnames(factors_M)]
    
    yhat<-function(x, i, j, Phe_M, factors_M) {
      train_envs = setdiff(colnames(factors_M), testEnvs)
      # 提取因子值和表型值
      x_train<-as.numeric(factors_M[i, train_envs])
      y_train<-as.numeric(Phe_M[j, train_envs])
      # 构建设计矩阵：列为 x^2, x, 1
      M_f<-cbind(x_train^2, x_train, 1)
      # 计算回归系数并预测
      coef<-ginv(M_f) %*% y_train
      y<-c(x^2, x, 1) %*% coef
      return(y)
    }
    
    y<-matrix(0, nrow = nrow(factors0), ncol = nrow(Phe0))
    y1<-matrix(0, nrow = nrow(factors0), ncol = nrow(Phe0))
    C1<-NULL
    for (i in 1:nrow(factors0)){
      for (j in 1:nrow(Phe0)){
        y[i,j]<-yhat(factors0[i, Envs], i, j, Phe_M, factors_M)
        y1[i,j]<-yhat(factors_M[i, testEnvs], i, j, Phe_M, factors_M)
      }
      C1[i]<-cor(y1[i,], Phe0[, testEnvs])
    }
    rownames(y) = rownames(factors0)
    colnames(y) = Phe0$line_code
    
    # Fixed Threshold
    yy<-data.frame(y[which(C1 >= C_cutoff), , drop = FALSE])
    yy<-yy[apply(yy, 1, function(row) all(row >= range_min & row <= range_max)), ]
    
    df_phe_PreEnv<-rbind(df_phe_PreEnv, yy)
  }
  
  Phe_PreEnv<-data.frame(line_code = colnames(df_phe_PreEnv),
                         Predicted = colMeans(df_phe_PreEnv), 
                         Traits = trait)
  return(Phe_PreEnv)
}
