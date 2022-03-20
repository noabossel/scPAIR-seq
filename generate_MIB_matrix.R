
mut_barcodes <- c("AGGCGGACAT","GACCGAACTG",
                  "TACCACACGG","CACTGGCGAA", "GAACTGTCCG","CATCGGCAGA",
                  "GGCGGATATG","CACCTACGTC", "GGTACCAGCA","TCGCACCGTA",
                  "GAGTTCGTCC","CTGGAGCTTC","ATGGACCAGG","ACGGTGTTCG",
                  "CAGCGCTCAT","CAATGGCACC","GCTTCACGAC","CAGAACGTCC",
                  "GCACTTGTCG","GCTGCGGATA","TCTTGGCGGA","CCGGTACCTT",
                  "ACGGTTCGCA","CGGACAGTTG","GGCCAGTCAA")

mut_strains <- c("WT","sspH2",
                 "sifB","sifA","spvB","gtgA",
                 "sseJ","srfJ","sseI","sseG",
                 "pipB2","gtgE","gogB","sopD2",
                 "sseL","steC","steB","sseF",
                 "steA","steE","sseK1","pipB",
                 "slrP","ssaV","spvD")

# load CEL-SEQ 96 barcodes
celseq_samples <- read.table("celseq_bc.txt")
names(celseq_samples) <- c("id","bc")
celseq_samples$name <- paste0("bc_", sprintf("%02d", 1:nrow(celseq_samples)))

# GFP sequence:
gfp_seq_1 <- "TGGGATTACACATGGCATGGATGAACTATACAAATAA" 
gfp_seq_2 <- "TGGGATCACACATGGCATGGATGAACTATACAAATAA" #for sseJ
GFP_mismatch <- 1

# load fastq files: 
fq_r1 <- readFastq("ex_R1.fastq")
fq_r2 <- readFastq("ex_R2.fastq")

#assign reads to CEL SEQ barcodes
cs_bc_in_fq <- substr(sread(fq_r1), start = 7, stop = 12)
celseq_counts <- c()
for(j in c(1:nrow(celseq_samples))){
  crnt_CS_reads <- (cs_bc_in_fq == celseq_samples[j,"bc"])
  assign(paste0("cs_",j,"_r1"), sread(fq_r1)[crnt_CS_reads])
  assign(paste0("cs_",j,"_r2"), sread(fq_r2)[crnt_CS_reads])
  #print(length(get(paste0("cs_",j,"_r1"))),quote = F)
  celseq_counts <- c(celseq_counts, length(get(paste0("cs_",j,"_r1"))))
}

used_cs_barcodes <- c(1:96)
all_96_cells_bc_1mm <- data.frame(matrix(0, ncol = length(mut_strains), 
                                         nrow = length(used_cs_barcodes)))

# iterate over all used CEL-SEQ barcodes
for(cs in used_cs_barcodes){
  
  # get CEL-SEQ to analyze
  crnt_smpl_r1 <- get(paste0("cs_",cs,"_r1"))
  crnt_smpl_r2 <- get(paste0("cs_",cs,"_r2"))
  
  # filter for GFP with specific number of mutations
  gfp_filter <- as.vector(t(adist(gfp_seq_1, substr(crnt_smpl_r2,1,37))))==GFP_mismatch | as.vector(t(adist(gfp_seq_2, substr(crnt_smpl_r2,1,37))))==GFP_mismatch
  crnt_smpl_r1_gfp <- crnt_smpl_r1[gfp_filter]
  crnt_smpl_r2_gfp <- crnt_smpl_r2[gfp_filter]
  
  r1_seqs <- as.character(crnt_smpl_r1_gfp)
  r2_seqs <- as.character(crnt_smpl_r2_gfp)
  
  #break R1 ad R2 reads to its components
  reads_df <- data.frame(UMI = substr(r1_seqs,1,6),
                         CEL_SEQ = substr(r1_seqs,7,12),
                         BARCODE = substr(r2_seqs,38,47))
  
 
  # check number of each mutant barcode in current sample 
  all_BC_df <- data.frame(CELSEQ = cs,
                          STRAIN = mut_strains,#all_strains,
                          BCS = mut_barcodes,#all_barcodes,
                          #PIE_TOT = length(r1_seqs),
                          stringsAsFactors = F)
  
  # find barcode with perfect match OR max of 1 mismatch
  all_BC_df$DETECTED_0 <- sapply(all_BC_df$BCS, function(x) sum(0 == adist(x,reads_df$BARCODE)))
  all_BC_df$DETECTED_1 <- sapply(all_BC_df$BCS, function(x) sum(1 >= adist(x,reads_df$BARCODE)))
  
  #add count of unknown barcodes as 'OTHER' 
  all_BC_df <- rbind(all_BC_df, c(cs, "OTHER","-",nrow(reads_df),nrow(reads_df)-sum(all_BC_df$DETECTED_0), nrow(reads_df)-sum(all_BC_df$DETECTED_1)))
  all_BC_df$DETECTED_0 <- as.numeric(all_BC_df$DETECTED_0)
  all_BC_df$DETECTED_1 <- as.numeric(all_BC_df$DETECTED_1)

  all_96_cells_bc_1mm[cs,] <- head(all_BC_df$DETECTED_1, -1)  
}

colnames(all_96_cells_bc_1mm) <- head(all_BC_df$STRAIN, -1)
write.table(all_96_cells_bc_1mm,
            "output_files/MIB_matrix.txt",
            quote = F, row.names = T) 



