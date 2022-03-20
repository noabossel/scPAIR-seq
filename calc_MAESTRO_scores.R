rm(list = ls())
library(data.table)
library("kdensity")
library(proxy)

source('MAESTRO_functions.R')

#######data for random#######
sce_all <- readRDS('sce_norm.RDS')

#separate cells for real (cells infected with single mutant) and random data (cells without detection of the invading bacteria used for the random model)
sce_for_random <- sce_all[,grepl('Removed|Negative|_Multi', sce_all$BACT_MUTANT)]
sce_single <- sce_all[,!grepl('Removed|Negative|_Multi', sce_all$BACT_MUTANT)]
n_single_cells <- ncol(sce_single) #757 cells infected with single mutant
mut_amount <- as.array((table(sce_single$BACT_MUTANT)))
mutant_vec <- sort(names(mut_amount), decreasing=T ) 
mutant_sort <- names(sort(mut_amount))

mut_col_names <- c()
for(mut in names(mut_amount)){
  mut_col_names <- c(mut_col_names, rep(mut, mut_amount[mut]))
}

#remove mito genes
table(grepl('^mt-', row.names(sce_for_random)))
rownames(sce_for_random)[grepl('^mt-', row.names(sce_for_random))]
mito_remove <- grep('^mt-', row.names(sce_for_random))

sce_for_random <- sce_for_random[-as.numeric(mito_remove),]

log_exp_for_random <-logcounts(sce_for_random)
lin_exp_for_random <- (2^log_exp_for_random)-1

n_random_cells <- ncol(sce_for_random) #2109
indexs_to_sample <- c(1:n_random_cells)

######real data#######
real_log_exp <- read.table('sce_log_counts_by_mut.txt',stringsAsFactors=F)
real_log_exp <- real_log_exp[,!grepl('.*_Multi', colnames(real_log_exp))]
colnames(real_log_exp) <- gsub('\\.\\d+','',colnames(real_log_exp))

real_lin_exp <- (2^real_log_exp)-1 #27853   757

######random sampling#########
set.seed(10)
indexs_df <- random_indexs_sampling(1000, indexs_to_sample, mut_amount, mut_col_names, "first")

##########CALCULATE MAESTRO SCORES#########

#for real
real_maestro_scores <- create_MAESTRO_scores(1.5, real_log_exp, real_lin_exp, log_mean_lin=F, mut_amount, cells_cutoff=0.3, low_umi_cutoff=2)
write.table(real_maestro_scores,"output_files/maestro_scores_real_data.txt", quote = F, sep = "\t")

#summary(real_maestro_scores$DISTANCE)

#generate random MAESTRO scores (null distribution per mutant)
rand_maestro_scores <- dist_hist_data_for_random_2(indexs_df, 1.5, F, log_exp_for_random, lin_exp_for_random, "output_files/random", mut_amount, cells_cutoff=0.3, low_umi_cutoff=2)

summary(rand_maestro_scores$DISTANCE)






