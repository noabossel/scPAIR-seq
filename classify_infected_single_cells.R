
#### Start from raw MIB matrix: ####
#==================================#
## choose parameters for detection analysis:
#crnt_data <- read.table("ex_MIB_matrix.txt", header = TRUE)
#rm_th <- 2^7

## remove cells with size under threshold
#rm_indices <- rowSums(crnt_data) < rm_th
#crnt_data <- crnt_data[!rm_indices,]
#print( dim(crnt_data) )

## create Seurat object with HTO data (our MIB):
#s_objct <- CreateSeuratObject(t(crnt_data),
#                              assay = 'counts')

#s_objct <- CreateSeuratObject(t(crnt_data),
#                              assay = 'HTO')

## normalize data by centered log ratio (geometric mean)
#s_objct <- NormalizeData(s_objct, assay = "HTO", normalization.method = "CLR")

## remove batch effect by plate
#plate_data <- readRDS("plate_data.RDS")
#satija_data <- t(s_objct@assays$HTO@data)
#satija_data_be <- t(removeBatchEffect(t(satija_data), plate_data[!rm_indices]))
#crnt_data <- 2^(satija_data_be)


#### Start from normalized MIB matrix: ####
#=========================================#
crnt_data <- read.table("ex_norm_MIB_matrix.txt", header = TRUE)
print( dim(crnt_data) )

#Perform Quantile Sweep to choose q that maximise singlets
bar.table <- crnt_data
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.01)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

## Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)

## Finalize round 1 classifications, remove negative cells
opt_q <- findQ(threshold.results1$res, threshold.results1$extrema)
round1.calls <- classifyCells(bar.table, q=opt_q)

write.table(round1.calls,
            "output_files/cells_detection.txt",
            quote = F, row.names = T) 





