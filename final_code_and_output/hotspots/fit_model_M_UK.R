library(data.table)
library(purrr)
library(parallel)
source("fns_sim.R")
#chr <- commandArgs(trailingOnly = TRUE)
set.seed(26)

# Detect the number of available cores
num_cores <- 12

read.tracts <- function(chr) {
  str <- paste0("/projects/browning/brwnlab/sharon/for_nobu/gc_length/ukbiobank/chr", chr, ".ibdclust2cM_trim1_combinedoffsets_v6.inf_obs_tracts2")
  # str <- paste0("050124//UK_biobank//chr", chr, ".ibdclust2cM_trim1_combinedoffsets_v6.inf_obs_tracts2")
  df <- read.table(str)
  return(df)
}

read.MAF <- function(chr) {
  str <- paste0("/projects/browning/brwnlab/sharon/for_nobu/gc_length/ukbiobank/chr", chr, ".allregions.pmaf.gz")
  # str <- paste0("050124//UK_biobank//chr", chr, ".allregions.pmaf.gz")
  df <- read.table(str)
  return(df)
}

tracts.df.0.5.lst <- lapply(1:22, read.tracts)
MAF.df.lst <- lapply(1:22, read.MAF)

boot_MLE_M <- function(x, l_lst, psi_lst, M, remove.singletons) {
  ind <- sample(1:length(l_lst), replace = TRUE)
  l_lst_boot <- l_lst[ind]
  psi_lst_boot <- psi_lst[ind]
  
  if (remove.singletons) { 
    optim.out.geom <- optim(0.005, calc_neg_log_lik, psi_lst = psi_lst_boot, pL = pL_geom_1M,
                            l_lst = l_lst_boot, M = M, lower = 0.0001, upper = 0.05, method = "Brent")}
  else { 
    optim.out.geom <- optim(0.005, calc_neg_log_lik, psi_lst = psi_lst_boot, pL = pL_geom_M,
                            l_lst = l_lst_boot, M = M, lower = 0.0001, upper = 0.05, method = "Brent")}
  
  print("finished 1 bootstrap iteration")
  return(list(optim.out.geom))
}

### For the UK Biobank data, we want to group all of the tract lengths and psi values from the 22 chromosomes. 
### So we define a function for this.
get_l_psi <- function(MAF.df, df, M, region, MAF.ceil, remove.singletons) {
  
  # to make indexing easier
  MAF.df$pos <- MAF.df$V2 + 1
  psi <- numeric(max(MAF.df$pos))
  psi[c(MAF.df$pos)] <- 2*MAF.df$V3*(1-MAF.df$V3)
  
  exc <- MAF.df[MAF.df$V3 < 0.05 | MAF.df$V3 > MAF.ceil ,]
  psi[exc$pos] <- 0
  # 
  # print("make sure this is 0: ")
  # print(psi[exc$pos])
  # 
  ### next, we calculate the psi for each tract using a region of 1000
  # get obs tracts
  df$V1 <- df$V1 + 1
  df$V2 <- df$V2 + 1
  # keep <- which( df$V2 - df$V1 + 1 <= M & df$V2 - df$V1 + 1 > 1 )
  if (remove.singletons) {keep <- which( df$V2 - df$V1 + 1 <= M & df$V2 - df$V1 + 1 > 1)}
  else {keep <- which( df$V2 - df$V1 + 1 <= M)}
  
  df <- df[keep, ]
  
  print(nrow(df))
  
  tracts_lst <- split(df, seq(nrow(df)))
  psi_lst <- lapply(tracts_lst, calculate_psi, psi = psi, region = region, length_chrom = length(psi), debias = FALSE)
  l_lst <- lapply(tracts_lst, function(x) {x[2]-x[1]+1})
  print("finished processing data for 1 chromosome")
  
  return(list(l_lst, psi_lst))
}

fit_model_M <- function(l_lst, psi_lst, M, remove.singletons) {
  
  ### maximum likelihood
  if (remove.singletons) {
    optim.out.geom <- optim(0.005, calc_neg_log_lik, psi_lst = psi_lst, pL = pL_geom_1M,
                            l_lst = l_lst, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  }
  else {
    optim.out.geom <- optim(0.005, calc_neg_log_lik, psi_lst = psi_lst, pL = pL_geom_M,
                            l_lst = l_lst, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  }
  
  res.boot <- mclapply(1:500, boot_MLE_M, l_lst, psi_lst, M, remove.singletons, mc.cores = num_cores)
  
  print("finished 1 seed")
  print(1/optim.out.geom$par)
  return(list(optim.out.geom, l_lst, psi_lst, res.boot))
}

l_psi_lst <- map2(MAF.df.lst, tracts.df.0.5.lst, get_l_psi, M = 1500, region = 5000, MAF.ceil = 0.5, remove.singletons = TRUE)
l_lst <- lapply(l_psi_lst, function(x) {x[[1]]}) %>% unlist()
psi_lst <- lapply(l_psi_lst, function(x) {x[[2]]}) %>% unlist()
print("number of obs. tracts: ")
print(length(l_lst))
print("length of psi: ")
print(length(psi_lst))

res.0.5.1M <- fit_model_M(l_lst, psi_lst, 1500, TRUE)
file_name <- paste0("res.UK_Biobank.1M.1500.region.5000.geom.ibdclust2cM.MAF.0.5.boot.grouped.keep.ends.rds")
saveRDS(res.0.5.1M, file_name)