library(data.table)
library(proxy)

#########general functions###############

#'@param- weight- numeric the weight for the mean
#'@param- log_exp- numeric matrix- log_norm exp
#'@param- lin_exp- numeric matrix- lin_norm exp
#'@param- log_mean_lin- boolean default=F if to use log(mean(linear)) or mean(log) for the MAESTRO score
#'@param- mut_amount- table- the mut names with cell amount
#'@param- cells_cutoff- numeric- default 0.3- for choosing the genes, to pass gene needed at least cells_cutoff*mut_cell_amount cells that express low_umi_cutoff or more UMIs 
#'@param- low_umi_cutoff- numeric- default 2- for choosing the genes, to pass gene needed at least cells_cutoff*mut_cell_amount cells that express low_umi_cutoff or more UMIs 
#'@return data frame of the MAESTRO scores
create_MAESTRO_scores = function(weight, log_exp, lin_exp, log_mean_lin=F, mut_amount, cells_cutoff=0.3, low_umi_cutoff=2){  
  cv2_df_list <- cv2_df_for_hist_2(log_exp, lin_exp, mut_amount, cells_cutoff, low_umi_cutoff) 
  cv2_with_WT_genes_df <- cv2_df_list[[1]]
  WT_with_muts_df_genes_df <- cv2_df_list[[2]]
  dist_hist_data <- create_distance_df(weight, cv2_with_WT_genes_df, WT_with_muts_df_genes_df, log_mean_lin)
  return(dist_hist_data)
}




#'@param- log_exp- numeric matrix- log_norm exp
#'@param- lin_exp- numeric matrix- lin_norm exp
#'@param- mut_amount- table- the mut names with cell amount
#'@param- cells_cutoff- numeric- default 0.3- for choosing the genes, to pass gene needed at least cells_cutoff*mut_cell_amount cells that express low_umi_cutoff or more UMIs 
#'@param- low_umi_cutoff- numeric- default 2- for choosing the genes, to pass gene needed at least cells_cutoff*mut_cell_amount cells that express low_umi_cutoff or more UMIs 
#'@return 

cv2_df_for_hist_2=function(log_exp, lin_exp, mut_amount, cells_cutoff=0.3, low_umi_cutoff=2){  
  WT_lin_exp <- lin_exp[, colnames(lin_exp) == "WT"]
  n_WT_cells <- mut_amount["WT"]
  #print(paste0("n_WT_cells=", n_WT_cells))
  WT_cells_cutoff<- cells_cutoff*n_WT_cells
  #print(paste0("WT_cells_cutoff=", WT_cells_cutoff))
  
  WT_lin_exp <- WT_lin_exp[apply(WT_lin_exp, 1, function(x) sum(x>=low_umi_cutoff)) >= WT_cells_cutoff, ]
  WT_pass_genes <- rownames(WT_lin_exp)
  #print(paste0("genes wt pass=", length(WT_pass_genes)))
  
  
  WT_log_exp <- log_exp[,colnames(log_exp)=="WT"]
  WT_log_exp <- WT_log_exp[WT_pass_genes, ]
  
  #get WT df
  WT_df <- get_cv2_df_per_mut("WT", WT_lin_exp, WT_log_exp, freq_threshold=2)
  
  #add genes from WT to MUT and genes from MUT to WT
  WT_with_muts_df_genes_df <- WT_df
  WT_with_muts_df_genes_df$is_log_modified <- FALSE
  WT_with_muts_df_genes_df$source <-"WT"
  
  cv2_with_WT_genes_df <- c() 
  
  for(i in c(2:length(mutant_vec))){
    crnt_ptrn <- mutant_vec[i]
    print(crnt_ptrn)
    crnt_lin_exp <- lin_exp[,colnames(lin_exp) == crnt_ptrn]
    n_mut_cells <- mut_amount[crnt_ptrn]
    mut_cells_cutoff<- cells_cutoff*n_mut_cells
    print(paste0("n_mut_cells=", n_mut_cells))
    print(paste0("mut_cells_cutoff=", mut_cells_cutoff))
    
    crnt_lin_exp <- crnt_lin_exp[apply(crnt_lin_exp, 1, function(x) sum(x>=low_umi_cutoff)) >= mut_cells_cutoff, ]
    mut_pass_genes <- rownames(crnt_lin_exp)
    print(paste0("genes mut pass=", length(mut_pass_genes)))
    
    cv2_mut_WT <- add_WT_genes_to_MUT_and_opp_to_cv2_df(crnt_ptrn, mut_pass_genes, WT_pass_genes, lin_exp, log_exp)
    cv2_mut_with_WT_genes <- cv2_mut_WT[[1]]
    cv2_WT_with_mut_genes <- cv2_mut_WT[[2]]
    
    #concate to dataframe  
    cv2_with_WT_genes_df <- rbind(cv2_with_WT_genes_df, cv2_mut_with_WT_genes)
    WT_with_muts_df_genes_df <- rbind(WT_with_muts_df_genes_df, cv2_WT_with_mut_genes)
  }
  #return(cv2_df)
  return(list(cv2_with_WT_genes_df, WT_with_muts_df_genes_df))
}



#'@param- weight- numeric the weight for the mean
#'@param- cv2_with_WT_genes_df- data frame- mean, cv2 df of the mut 
#'@param- WT_with_muts_df_genes_df- data frame- mean, cv2 df of the WT 
#'@param- log_mean_lin- boolean default=F if to use log(mean(linear)) or mean(log) for the MAESTRO scores
#'@return data frame of the MAESTRO scores

create_distance_df = function(weight, cv2_with_WT_genes_df, WT_with_muts_df_genes_df, log_mean_lin=F){
  dist_hist_data <- data.frame()
  
  for(i in c(1:length(mutant_vec))){
    crnt_mut <- mutant_vec[i]
    if(crnt_mut == 'WT'){next}
    
    #get data of current mutant
    crnt_mut_df <- cv2_with_WT_genes_df[cv2_with_WT_genes_df$MUTANT==crnt_mut,]
    WT_with_crnt_mut_df <- WT_with_muts_df_genes_df[WT_with_muts_df_genes_df$source %in% c(crnt_mut, "WT"), ]
    rownames(crnt_mut_df) <- crnt_mut_df$GENE
    rownames(WT_with_crnt_mut_df) <- WT_with_crnt_mut_df$GENE
    WT_with_crnt_mut_df <- WT_with_crnt_mut_df[rownames(crnt_mut_df), ]
    
    #calculate distance between points of mutant and WT
    if(log_mean_lin){
      log2_crnt_mut_df <- log2(crnt_mut_df[, c('ptrn_means','ptrn_cv2')])
      log2_crnt_mut_df$weight_ptrn_means <- weight*(log2_crnt_mut_df$ptrn_means)
      log2_WT_with_crnt_mut_df <- log2(WT_with_crnt_mut_df[, c('ptrn_means','ptrn_cv2')])
      log2_WT_with_crnt_mut_df$weight_ptrn_means <- weight*(log2_WT_with_crnt_mut_df$ptrn_means)
      crnt_mut_dist <-  dist(log2_crnt_mut_df[,c('weight_ptrn_means','ptrn_cv2')],
                             log2_WT_with_crnt_mut_df[,c('weight_ptrn_means','ptrn_cv2')],
                             method = "euclidean", pairwise = TRUE)
    }
    
    
    #log2_crnt_mut_df <- log2(crnt_mut_df[, c('ptrn_means','ptrn_cv2')])
    else{
      crnt_mut_df$weight_ptrn_log_means <- weight*(crnt_mut_df$ptrn_log_means)
      WT_with_crnt_mut_df$weight_ptrn_log_means <- weight*(WT_with_crnt_mut_df$ptrn_log_means)
      
      crnt_mut_dist <-  dist(crnt_mut_df[,c('weight_ptrn_log_means','ptrn_log_cv2')],
                             WT_with_crnt_mut_df[,c('weight_ptrn_log_means','ptrn_log_cv2')],
                             method = "euclidean", pairwise = TRUE)
    }
    #convert class
    crnt_mut_dist <- as.numeric(crnt_mut_dist)
    
    #keep all relevant data
    crnt_hist_data <- data.frame(MUTANT_REF = 'WT',
                                 MUTANT = crnt_mut,
                                 GENE = rownames(WT_with_crnt_mut_df),
                                 REF_MOD = WT_with_crnt_mut_df$is_log_modified,
                                 MUTANT_MOD = crnt_mut_df$is_log_modified,
                                 DISTANCE = crnt_mut_dist,
                                 crnt_mut_df[,c('ptrn_log_means','ptrn_log_cv2','ptrn_exp_freq', 'ptrn_means', 'ptrn_cv2')],
                                 #P99 = unname(quantile(crnt_mut_dist,0.88)),
                                 REF_MEAN = WT_with_crnt_mut_df$ptrn_log_means,
                                 #REF_weight_MEAN = WT_with_crnt_mut_df$weight_ptrn_means,
                                 REF_CV2 = WT_with_crnt_mut_df$ptrn_log_cv2,
                                 REF_MEAN_lin = WT_with_crnt_mut_df$ptrn_means,
                                 REF_CV2_lin = WT_with_crnt_mut_df$ptrn_cv2,
                                 REF_EXP_FREQ = WT_with_crnt_mut_df$ptrn_exp_freq)#,
    #source =crnt_mut_dist$source)
    
    dist_hist_data <- rbind(dist_hist_data, crnt_hist_data)
  }
  return(dist_hist_data)
}



#'create data frame of the mean, cv2 for MUTANT
#'@param- mut- string-  MUTANT name
#'@param- mut_lin_exp- numeric matrix- lin_norm exp
#'@param- mut_log_exp- numeric matrix- log_norm exp
#'@param- freq_threshold - numeric- (need to be equal to low_UMI from the other functions)
#'@return data frame of the mean, cv2 for MUTANT
get_cv2_df_per_mut = function(mut, mut_lin_exp, mut_log_exp, freq_threshold=2){
  #get mutant mean and variance
  mut_means <- rowMeans(mut_lin_exp)
  mut_vars <- rowVars(as.matrix(mut_lin_exp))
  mut_log_means <- rowMeans(mut_log_exp)
  mut_log_vars <- rowVars(as.matrix(mut_log_exp))
  
  #modify means to avoid NA values
  is_log_modified <- mut_log_means==0
  print(table(is_log_modified))
  is_modified <- mut_means==0
  #print(table(is_modified))
  #set 0 to minimum value of current mutant
  mut_vars[mut_vars==0] <- min(mut_vars[mut_vars!=0])
  mut_log_vars[mut_log_vars==0] <- min(mut_log_vars[mut_log_vars!=0])
  
  mut_means[mut_means==0] <- min(mut_means[mut_means!=0])
  mut_log_means[mut_log_means==0] <- min(mut_log_means[mut_log_means!=0])
  
  
  # calc coefficient of variance
  mut_cv2 <- sqrt(mut_vars)/mut_means
  mut_log_cv2 <- sqrt(mut_log_vars)/mut_log_means 
  
  
  #freq of expressing cells
  mut_exp_freq <- apply(mut_lin_exp, 1, function(x) sum(x>=freq_threshold)/length(x))
  
  
  #keep data
  mut_df <- data.frame(ptrn_means=mut_means,
                       ptrn_log_means=mut_log_means,
                       ptrn_vars=mut_vars,
                       ptrn_log_vars=mut_log_vars,
                       ptrn_cv2=mut_cv2,
                       ptrn_log_cv2=mut_log_cv2,
                       ptrn_exp_freq=mut_exp_freq,
                       MUTANT=mut,
                       is_log_modified = is_log_modified,
                       GENE=rownames(mut_lin_exp))
  return(mut_df)
}




#'@param- mut- string- MUTANT name
#'@param- mut_genes- string vector- MUTANT genes that passed
#'@param- WT_genes- string vector- WT genes that passed
#'@param- log_exp- numeric matrix- log_norm exp
#'@param- lin_exp- numeric matrix- lin_norm exp
#'@return data frame of the MAESTRO scores
add_WT_genes_to_MUT_and_opp_to_cv2_df = function(mut, mut_genes, WT_genes, lin_exp, log_exp){
  only_mut_genes <- setdiff(mut_genes, WT_genes)
  only_WT_genes <- setdiff(WT_genes, mut_genes)
  #print(paste0("intersect genes=", length(intersect(mut_genes, WT_genes))))
  #print(paste0("only ", mut, "genes=", length(only_mut_genes)))
  #print(paste0("only WT genes=", length(only_WT_genes)))
  
  mut_lin_exp <- lin_exp[c(as.vector(mut_genes), as.vector(only_WT_genes)), colnames(lin_exp)==mut]
  mut_log_exp <- log_exp[c(as.vector(mut_genes), as.vector(only_WT_genes)), colnames(log_exp)==mut]
  WT_mut_lin_exp <- lin_exp[as.vector(only_mut_genes), colnames(lin_exp)=="WT"]
  WT_mut_log_exp <- log_exp[as.vector(only_mut_genes), colnames(log_exp)=="WT"]
  
  mut_WT_df <- get_cv2_df_per_mut(mut, mut_lin_exp, mut_log_exp, freq_threshold=2)
  WT_MUT_df <- get_cv2_df_per_mut("WT", WT_mut_lin_exp, WT_mut_log_exp, freq_threshold=2)
  
  source <- c(rep(mut, length(mut_genes)), rep("WT", length(only_WT_genes))) 
  WT_source <- c(rep(mut, length(only_mut_genes))) 
  
  mut_WT_df$source <- source
  WT_MUT_df$source <-WT_source
  
  print(dim(WT_MUT_df))
  
  return(list(mut_WT_df, WT_MUT_df))
  
}


#'save RDS and csv of the dataframe (each column is ona random sample, and one column of the mut "names")
#'@param- times- numeric- how many times to sample
#'@param- indexs_to_sample- numeric vector- sample from this vector
#'@param- mut_amount- table- the mut names with cell amount
#'@param- mut_col_names- string vector- mut names each mut*mut cell number
#'@param- dir- string- the directory for saving the files
#'@return data frame with diminution(length(mut_col_names), times+1), each column is the indexes of one sample round, and one column(MUTANT) is the "mut" 
random_indexs_sampling = function(times, indexs_to_sample, mut_amount, mut_col_names, dir){
  indexs_df <- data.frame(MUTANT=mut_col_names)
  for(i in c(1:times)){
    mut_indexs_vec <- c()
    remaining_indexs_to_sample <- indexs_to_sample
    for(mut in names(mut_amount)){
      mut_indexs <- sample(remaining_indexs_to_sample, size=mut_amount[mut], replace =F)
      remaining_indexs_to_sample <- setdiff(remaining_indexs_to_sample, mut_indexs)
      mut_indexs_vec <- c(mut_indexs_vec, mut_indexs)
    }
    indexs_df[, as.character(i)] <- mut_indexs_vec
  }
  return(indexs_df)
}


#'save RDS of MAESTRO data frames list 
#'@param- indexs_df- dataframe of the random indexes sampling- each column is one round of sampling, and the last column is the mutant "names"
#'@param- weight- numeric the weight for the mean
#'@param- log_mean_lin- boolean default=F if to use log(mean(linear)) or mean(log) for the MAESTRO scores
#'@param- log_exp_for_random- numeric matrix- the log_norm of the random cells
#'@param- lin_exp_for_random- numeric matrix- the lin_norm of the random cells
#'@param- dir- string- the directory for saving the files
#'@param- mut_amount- table- the mut names with cell amount
#'@param- cells_cutoff- numeric- default 0.3- for choosing the genes, to pass gene needed at least cells_cutoff*mut_cell_amount cells that express low_umi_cutoff or more UMIs 
#'@param- low_umi_cutoff- numeric- default 2- for choosing the genes as explained above 
#'@return data frame of the MUTANT, GENE, and MAESTRO scores (of all the 1000 together)
dist_hist_data_for_random_2=function(indexs_df, weight, log_mean_lin, log_exp_for_random, lin_exp_for_random, dir, mut_amount, cells_cutoff=0.3, low_umi_cutoff=2){
  dir.create(dir, recursive = T)
  dist_hists <- list()
  dist_hists_relevant <- list()
  for(i in as.character(c(1:(ncol(indexs_df)-1)))){
    print(paste0("starting round ", i))
    log_exp <- log_exp_for_random[, indexs_df[,i]]
    lin_exp <- lin_exp_for_random[, indexs_df[,i]]
    colnames(log_exp) <- mut_col_names
    colnames(lin_exp) <- mut_col_names
    dist_hist_data <- create_MAESTRO_scores(weight, log_exp, lin_exp, log_mean_lin, mut_amount, cells_cutoff, low_umi_cutoff)
    dist_hists[[i]] <- dist_hist_data
    dist_hists_relevant[[i]] <- dist_hist_data[, c("MUTANT", "GENE", "DISTANCE")]
  }
  
  all_dist_hist_relevant <- rbindlist(dist_hists_relevant)
  #saveRDS(all_dist_hist_relevant, paste0(dir, "/all_dist_hist_relevant.RDS"))
  #saveRDS(dist_hists,  paste0(dir, "/dist_hists_list.RDS"))
  return(all_dist_hist_relevant)
}
