# The script was used to obtain the optimal weights and select similar individuals for a given target.
# Notes:
# N: pool sizes
# predicted_train0: the predicted values of the traits in the training set  
# predicted_test: the predicted values of the traits in the test set  
# observeded_train: the observeded values of the traits in the training set
# observeded_test: the observeded values of the traits in the test set
# W: the optimal weights
# Pre_test: the normalized predicted values of the traits in the test set   
# Obs_test: the normalized observeded values of the traits in the test set
# n: pool size
# m: the number of random permutations
# target: a given target

prePro<-function(data){
 Mean<-apply(data, 2, mean)
 SD<-apply(data, 2, sd)
 value<-apply(data, 2, function(x){(x - mean(x)) / sd(x)})
 model = list(Mean = Mean, SD = SD, value = value)
 return(model)
}
pro<-function(w, a_M){
  pi_M<-NULL
  proi<-matrix(0, nrow = dim(a_M)[1], ncol = dim(a_M)[1])
  sf<-matrix(0, nrow = dim(a_M)[1], ncol = dim(a_M)[2])
  for(i in 1:dim(a_M)[1]){
    pi_M<-apply(a_M[, , i], 1, function(x){w %*% x})
    proi[i,]<-exp(pi_M) / (sum(exp(pi_M)))
    sf[i, ]<-proi[i, ]%*%a_M[, , i]
  }
  model = list(proi = proi, sf = sf)
  return(model)
}

ide_rate<-function(w, Obs_test, Pre_test, n, m){
  rate<-NULL
  for(k in 1:m){
    set.seed(k)
    ord<-sample(1:nrow(Obs_test))
    Pre_test_M<-Pre_test[ord, ]
    Obs_test_M<-Obs_test[ord, ]
    p<-array(0, dim = c(n, n, dim(Pre_test_M)[1] / n))
    for (j in 1:(dim(Pre_test_M)[1] / n)){
      pre_test1<-as.matrix(Pre_test_M[(n * (j - 1) + 1):(n * j), ])
      obs_test1<-as.matrix(Obs_test_M[(n * (j - 1) + 1):(n * j), ])
      Dis<-abs(obs_test1 - pre_test1)
      m_M<-matrix(1, nrow = dim(obs_test1)[1], ncol = 1)
      a_M<-array(0, dim = c(dim(pre_test1)[1], dim(pre_test1)[2], dim(pre_test1)[1]))
      for(i in 1:dim(pre_test1)[1]){
        a_M[, , i]<-abs(m_M %*% as.numeric(pre_test1[i, ]) - obs_test1)
      }
      p[, , j]<-pro(w, a_M)$proi
    }
    a<-NULL
    for(j in 1:dim(Pre_test_M)[1] / n){
      a[j]<-sum(apply(p[, , j], 1, which.max) == c(1:n))
    }
    rate[k]<-sum(a) / (floor(dim(Pre_test_M)[1] / n) * n)
  }
  return(mean(rate))
}

###### Learn optimal weights
Weight_res<-function(predicted_train0, observeded_train, names_trait, b){
  
  Ib<-which(predicted_train0[nrow(predicted_train0), ] >= b)
  names_trait_M<-names_trait[Ib]
  predicted_train<-predicted_train0[1:nrow(observeded_train), ]
  
  pre_Pro<-prePro(predicted_train[, Ib])
  obs_Pro<-prePro(observeded_train[, Ib])
  
  predicted_train_M<-pre_Pro$value
  observeded_train_M<-obs_Pro$value
  
  Pre_train<-as.matrix(predicted_train_M)
  Obs_train<-as.matrix(observeded_train_M)
  Dis<-abs(Pre_train - Obs_train)
  
  m_M<-matrix(1, nrow = dim(Obs_train)[1], ncol = 1)
  a_M<-array(0, dim = c(dim(Pre_train)[1], dim(Pre_train)[2], dim(Pre_train)[1]))
  
  
  for(i in 1:dim(Pre_train)[1]){
    a_M[, , i]<-abs(m_M %*% as.numeric(Pre_train[i, ]) - Obs_train)
  }
  
  fr<-function(w){
    pi_M<-NULL
    proi<-matrix(0, nrow = dim(a_M)[1], ncol = dim(a_M)[1])
    sf<-matrix(0, nrow = dim(a_M)[1], ncol = dim(a_M)[2])
    for(i in 1:dim(a_M)[1]){
      pi_M<-apply(a_M[, , i], 1, function(x){w %*% x})
      proi[i, ]<-exp(pi_M) / (sum(exp(pi_M)))
      sf[i, ]<-proi[i, ] %*% a_M[, , i]
    }
    f<-sum(log(diag(proi)))
    return(f)
  }
  gr<-function(w){
    pi_M<-NULL
    proi<-matrix(0, nrow = dim(a_M)[1], ncol = dim(a_M)[1])
    sf<-matrix(0, nrow = dim(a_M)[1], ncol = dim(a_M)[2])
    for(i in 1:dim(a_M)[1]){
      pi_M<-apply(a_M[, , i], 1, function(x){w %*% x})
      proi[i, ]<-exp(pi_M) / (sum(exp(pi_M)))
      sf[i, ]<-proi[i, ] %*% a_M[, , i]
    }
    g<-apply(Dis, 2, sum) - apply(sf, 2, sum)
    return(g)
  }
  res<-constrOptim(rep(0, ncol(Obs_train)), fr, gr, method = "BFGS", control = list(fnscale = -1, maxit = 1500), ui = rbind(rep(1, ncol(Obs_train))), ci = rep(-20, 1))
  w<-res$par
  W<-(-w)
  W_matrix<-data.frame(trait = names_trait_M, weight = W)
  model=list(W_matrix = W_matrix, pre_mean = pre_Pro$Mean, pre_sd = pre_Pro$SD, obs_mean = obs_Pro$Mean, obs_sd = obs_Pro$SD)
  return(model)
}

###### Test the top model
Test_top_acc<-function(predicted_test, observeded_test, W, names_trait, N, m){ 
  Ib<-match(W[, 1], names_trait)
  predicted_test_M<-prePro(predicted_test[, Ib])$value
  observeded_test_M<-prePro(observeded_test[, Ib])$value
  
  w<-(-W[, 2])
  Ide_rate<-NULL
  for (i in 1:length(N)){
    Ide_rate[i]<-ide_rate(w, observeded_test_M, predicted_test_M, N[i], m)
  }
  return(Ide_rate)
}

Demo_P<-function(w, Obs_test, Pre_test, n, m){
  for(k in 1:m){
    set.seed(123)
    ord<-sample(1:nrow(Obs_test))
    Pre_test_M<-Pre_test[ord, ]
    Obs_test_M<-Obs_test[ord, ]
    p<-array(0, dim = c(n, n, dim(Pre_test_M)[1] / n))
    for (j in 1:(dim(Pre_test_M)[1] / n)){
      pre_test1<-as.matrix(Pre_test_M[(n * (j - 1) + 1):(n * j), ])
      obs_test1<-as.matrix(Obs_test_M[(n * (j - 1) + 1):(n * j), ])
      Dis<-abs(obs_test1 - pre_test1)
      m_M<-matrix(1, nrow = dim(obs_test1)[1], ncol = 1)
      a_M<-array(0, dim = c(dim(pre_test1)[1], dim(pre_test1)[2], dim(pre_test1)[1]))
      for(i in 1:dim(pre_test1)[1]){
        a_M[, , i]<-abs(m_M %*% as.numeric(pre_test1[i, ]) - obs_test1)
      }
      p[, , j]<-pro(w, a_M)$proi
    }
  }
  similarity_matrix<-p[, , 1]
  rownames(similarity_matrix) = colnames(similarity_matrix) = paste0('id', 1:5)
  return(similarity_matrix)
}

###### Select similar individuals for a given target
Tpro<-function(w, ta_M){
	pi_M<-NULL
	tproi<-matrix(0, nrow = 1, ncol = dim(ta_M)[1])
	tpi_M<-apply(ta_M, 1, function(x){w%*%x})
	tproi<-exp(tpi_M) / sum(exp(tpi_M))
  return(tproi)
}
Top_target<-function(target0, predicted_test, names_test, names_trait, pre_train_mean, pre_train_sd, obs_train_mean, obs_train_sd, selection_ratio, improve_ratio, improve_trait, W){
  Ib<-match(W[, 1], names_trait)
  target1<-target0
  if (!("NA" %in% improve_trait) || length(improve_trait) > 1){
    target1[improve_trait]<-(1 + improve_ratio)*target0[improve_trait]
    target<-target1[Ib]
  }else{
    target<-target1[Ib]
  }
  w<-(-W[,2])
  Pre_test<-predicted_test[, Ib, drop = F]
  target_M<-(target-pre_train_mean) / pre_train_sd
  predicted_test_M<-(Pre_test - matrix(1, nrow = nrow(Pre_test), ncol = 1) %*% pre_train_mean) / matrix(1, nrow = nrow(Pre_test), ncol = 1) %*% pre_train_sd
  Target<-as.matrix(rep(1, nrow(predicted_test_M))) %*% t(as.matrix(target_M))
  ta_M<-abs(Target - predicted_test_M)
  similarity0<-Tpro(w, ta_M)
  select_order<-match(sort(similarity0, decreasing = TRUE), similarity0)
  select<-select_order[1:(nrow(Pre_test) * selection_ratio)]
  names_select<-names_test[select]
  values_select<-Pre_test[select, ]
  model=list(names_select = names_select, values_select = values_select, similarity0 = similarity0)
  return(model)
}

Gpro<-function(Ind_S, G0, index_test, index_train){
  G<-apply(G0[Ind_S, ], 2, sum)
  G_M<-G
  G_M<-as.matrix(G_M)
  G_test<-as.matrix(G_M[index_test, ])
  G_train<-as.matrix(G_M[index_train, ])
  model = list(G_test = G_test, G_train = G_train)
  return(model)
}

Amatrix<-function(trait_train, Index_g){
  trait_trainN0<-trait_train[, -1]
  trait_trainN<-apply(trait_trainN0, 2, as.numeric)
  a<-apply(trait_trainN[, Index_g], 2, function(x){(x-mean(x)) / sd(x)})
  A_M<-a%*%t(a)
  A_N<-A_M
  for (i in 1:nrow(A_M)){
    for (j in 1:nrow(A_M)){
      A_N[i, j]<-A_M[i, j] / (sqrt(A_M[i, i]) * sqrt(A_M[j, j]))
    }
  }
  rownames(A_N)<-trait_train[[1]]
  colnames(A_N)<-trait_train[[1]]
  return(A_N)
}

###### select similar individuals for a given target
choose_material<-function(trait, trait_related, selection_ratio, improve_ratio, targetID){
  # Training set phenotype
  dis0<-Phe_Envs_Train %>%
    pivot_wider(names_from = Traits, values_from = Predicted)
  # Training set pseudo-true phenotype
  res_M0<-Phe_Envs_Train_Pseudo %>%
    pivot_wider(names_from = Traits, values_from = Predicted)
  # Test set phenotype
  Tdis0<-Phe_Envs_Test %>%
    pivot_wider(names_from = Traits, values_from = Predicted)
  Tres_M0<-Tdis0
  
  # 10-Fold Order
  set.seed(923)
  list_ID = as.character(dis0$line_code)
  list_ID = sample(list_ID)
  dis1<-dis0[match(list_ID, dis0$line_code), ]
  res_M1<-res_M0[match(list_ID, res_M0$line_code), ]
  # Cor
  Cor<-NULL
  Cor[1]<-"Cor"
  for (i in 2:ncol(res_M1)){
    Cor[i]<-cor(dis1[, i], res_M1[, i])
  }
  res_M1<-rbind(res_M1, Cor)
  
  # Extract chr and pos as matrices
  G_P<-rownames(Genotype)
  Ind_all<-str_match(G_P, "^chr(\\d+)\\.s_(\\d+)$")[, 2:3] %>% apply(2, as.numeric)
  # Extract index
  index_train<-match(dis1[[1]], colnames(Genotype))
  index_test<-match(Tdis0[[1]], colnames(Genotype))
  # Specifying trait phenotypes
  trait_test<-Tdis0
  trait_train<-dis1
  names_test<-Tdis0[, 1]
  
  dis_1<-apply(dis1[, 2:ncol(dis1)], 2, as.numeric)
  res_1<-apply(res_M1[, 2:ncol(res_M1)], 2, as.numeric)
  Tdis_1<-apply(Tdis0[, 2:ncol(Tdis0)], 2, as.numeric)
  Tres_1<-apply(Tres_M0[, 2:ncol(Tres_M0)], 2, as.numeric)
  
  Index_g2 = 1:(ncol(dis0) - 1)
  # Find the SNP site index that falls within the gene interval (S_pos)
  Gene_pos<-Gene_SNP_pos[, c('Chr', 'Start', 'End', 'GeneID')]
  Gene_pos<-Gene_pos[!duplicated(Gene_pos), ]
  
  S_pos<-integer()
  for (i in 1:nrow(Gene_pos)) {
    chr_i = Gene_pos$Chr[i]
    start_i = Gene_pos$Start[i]
    end_i = Gene_pos$End[i]
    
    S_pos_mid = which(Ind_all[, 1] == chr_i & Ind_all[, 2] >= start_i & Ind_all[, 2] <= end_i)
    S_pos = c(S_pos, S_pos_mid)
  }
  S_pos = unique(S_pos)
  
  # Extract the genotype matrix of the training set and the test set
  Geno_M<-t(Genotype[, 1:ncol(Genotype)])
  train_Geno<-Geno_M[index_train, ]
  test_Geno<-Geno_M[index_test, ]
  

  check_G0<-test_Geno[rownames(test_Geno) %in% targetID, S_pos, drop = FALSE]
  # Increase trait value
  S_pos_G<-NULL
  pv_G<-NULL
  
  improve_trait<-which(colnames(dis0) %in% trait)
  for (i in 1:length(S_pos)){
    Index_check_G<-which(train_Geno[, S_pos[i]] == check_G0[i])
    pv_G[i]<-t.test(dis1[Index_check_G, improve_trait], dis1[-Index_check_G, improve_trait], "greater")$p.value
    if (pv_G[i] < 0.05){
      S_pos_G<-c(S_pos_G, S_pos[i])
    }
  }
  
  
  list_SNP_ID = rownames(Genotype)[S_pos_G]
  SNP_pos<-data.frame(do.call(rbind, strsplit(sub("chr", "", list_SNP_ID), "\\.s_")))
  colnames(SNP_pos) = c('Chr', 'SNP_Pos')
  df_SNP0<-merge(SNP_pos, Gene_SNP_pos, by = c('Chr', 'SNP_Pos'), all.x = T)
  df_SNP<-df_SNP0[, c("Chr", "Start", "End", "GeneID")] %>%
    mutate(across(1:3, as.numeric)) %>%
    distinct()
  
  df_SNP0$SNP_ID = paste0("chr", df_SNP0$Chr, ".s_", df_SNP0$SNP_Pos)
  Genotype_New<-Genotype[rownames(Genotype) %in% df_SNP0$SNP_ID, ]
  
  Index_g2 = 1:(ncol(dis1) - 1)
  for (GeneNum in 1:nrow(df_SNP)) {
    current_gene = df_SNP[GeneNum, , drop = FALSE] 
    
    S_pos<-NULL
    for (i in 1:nrow(current_gene)) {
      Ind_chr = current_gene[i, 1]
      Ind_chos0<-Ind_all[Ind_all[, 1] == Ind_chr, ]
      
      if (Ind_chr > 1) {
        offset<-max(which(Ind_all[, 1] == (Ind_chr - 1)))
        S_pos<-c(S_pos, which((Ind_chos0[, 2] >= current_gene[i, 2]) & 
                                  (Ind_chos0[, 2] <= current_gene[i, 3])) + offset)
      } else {
        S_pos<-c(S_pos, which((Ind_chos0[, 2] >= current_gene[i, 2]) & 
                                  (Ind_chos0[, 2] <= current_gene[i, 3])))
      }
    }
    
    S_pos<-which(Genotype_New[, 1] %in% Genotype[S_pos, 1])
    G1_M<-Gpro(S_pos, Genotype_New, index_test, index_train)
    G_test<-G1_M$G_test
    G_train<-G1_M$G_train
    
    Am<-Amatrix(trait_train, Index_g2)
    num_CV_test<-round(nrow(G_train) / 10)
    A_train<-Am
    G_prediction_train<-matrix(0, nrow = nrow(G_train) + 11, ncol = ncol(G_train))
    for (i in 1:ncol(G_train)){
      G_prediction_train[, i]<-Pre_CV10(G_train[, i], A_train, num_CV_test)
    }
    G_train<-apply(G_train, 2, as.numeric)
    G_test<-apply(G_test, 2, as.numeric)
    
    gene_name = current_gene[1, 4]
    prediction_train_1<-cbind(res_1[, Index_g2], G_prediction_train[c(1:(nrow(res_1)-1), nrow(G_prediction_train)), 1])
    colnames(prediction_train_1)[ncol(prediction_train_1)] = gene_name
    
    trait_train_1<-cbind(dis_1[, Index_g2], G_train[, 1])
    colnames(trait_train_1)[ncol(trait_train_1)] = gene_name
    
    prediction_test_1<-cbind(Tdis_1[, Index_g2], G_test[, 1])
    colnames(prediction_test_1)[ncol(prediction_test_1)] = gene_name
    
    trait_test_1<-cbind(Tdis_1[, Index_g2], G_test[, 1])
    colnames(trait_test_1)[ncol(trait_test_1)] = gene_name
    
    if (GeneNum == 1) {
      df_prediction_train<-prediction_train_1
      df_trait_train<-trait_train_1
      df_prediction_test<-prediction_test_1
      df_trait_test<-trait_test_1
    }else{
      df_prediction_train<-cbind(df_prediction_train, prediction_train_1[, ncol(dis1), drop = FALSE])
      df_trait_train<-cbind(df_trait_train, trait_train_1[, ncol(dis1), drop = FALSE])
      df_prediction_test<-cbind(df_prediction_test, prediction_test_1[, ncol(dis1), drop = FALSE])
      df_trait_test<-cbind(df_trait_test, trait_test_1[, ncol(dis1), drop = FALSE])
    }
  }
  
  
  names_trait_1 = colnames(df_trait_train)
  
  m<-50
  Optimal_Weight_1<-Weight_res(df_prediction_train,
                               df_trait_train,
                               names_trait_1,
                               b = 0.2)
  save(Optimal_Weight_1, file = paste0("01_Weight_", trait, "_Late_NoSTD.RData"))
  
  Weight_1<-Optimal_Weight_1$W_matrix
  
  
  target0<-df_trait_test[Tdis0$line_code == targetID, Weight_1$trait]
  predicted_test<-df_prediction_test[Tdis0$line_code != targetID, Weight_1$trait, drop = FALSE]
  names_test<-Tdis0[Tdis0$line_code != targetID, 1][[1]]
  names_trait<-Weight_1$trait
  pre_train_mean<-Optimal_Weight_1$pre_mean
  pre_train_sd<-Optimal_Weight_1$pre_sd
  obs_train_mean<-Optimal_Weight_1$obs_mean
  obs_train_sd<-Optimal_Weight_1$obs_sd

  improve_trait<-which(colnames(predicted_test) %in% c(trait, trait_related))
  
  W<-Weight_1
  T_early<-Top_target(target0, predicted_test, names_test, names_trait, pre_train_mean, 
                      pre_train_sd, obs_train_mean, obs_train_sd, selection_ratio, improve_ratio, improve_trait, W)
  index_ordered<-order(T_early$values_select[, which(colnames(predicted_test) == trait)])
  index_name0<-T_early$names_select[index_ordered]
  T_name<-index_name0[1:100]
  Top50<-Tdis0[Tdis0$line_code %in% T_name, ] %>% 
    dplyr::arrange(desc(EW)) %>% 
    head(50)
  
  Top_ID = Top50$line_code
  return(Top_ID)
}

