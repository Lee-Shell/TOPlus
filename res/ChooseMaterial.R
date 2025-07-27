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
