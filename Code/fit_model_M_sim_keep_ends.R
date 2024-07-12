library(data.table)
library(purrr)
library(parallel)
source("fns_sim.R")
set.seed(26)

# Detect the number of available cores
num_cores <- 20

read.tracts <- function(seed, MAF) {
  if (MAF == 0.5) {
    str = paste0("/projects/browning/brwnlab/masakin/gene_conversion/sim5_data/sim5_seed", seed,
                 "_10Mb_n125000_ibdclust2cM_trim1_maf0.1_gcmaf0.05_9kbgaps_combinedoffsets_err0.0002_del1_1.0_421.inf_obs_tracts")
  }
  if (MAF == 0.4) {
    str = paste0("/projects/browning/brwnlab/masakin/gene_conversion/sim5_data/sim5_seed", seed,
                 "_10Mb_n125000_ibdclust2cM_trim1_maf0.1_gcmaf0.05_9kbgaps_combinedoffsets_err0.0002_del1_0.4_421.inf_obs_tracts")
  }
  if (MAF == 0.3) {
    str = paste0("/projects/browning/brwnlab/masakin/gene_conversion/sim5_data/sim5_seed", seed,
                 "_10Mb_n125000_ibdclust2cM_trim1_maf0.1_gcmaf0.05_9kbgaps_combinedoffsets_err0.0002_del1_0.3_421.inf_obs_tracts")
  }
  df <- read.table(str)
  return(df)
}

read.MAF <- function(seed) {
  str = paste0("/projects/browning/brwnlab/masakin/gene_conversion/sim5_data/sim5_seed", seed,
               "_10Mb_n125000.gtstats")
  df <- read.table(str)
  return(df)
}

tracts.df.list.0.5 <- lapply(1:20, read.tracts, MAF = 0.5)

MAF.df.list <- lapply(1:20, read.MAF)

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

fit_model_M <- function(MAF.df, df, M, region, MAF.ceil, remove.singletons) {
  
  # to make indexing easier
  MAF.df$pos <- MAF.df$V2 + 1
  psi <- numeric(max(MAF.df$pos))
  psi[c(MAF.df$pos)] <- 2*MAF.df$V11*(1-MAF.df$V11)
  
  exc <- MAF.df[MAF.df$V11 < 0.05 | MAF.df$V11 > MAF.ceil ,]
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
  
  print(l_lst %>% unlist())
  
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
  return(list(optim.out.geom, tracts_lst, psi_lst, res.boot))
}

res.0.5.1M <- pmap(list(MAF.df.list, tracts.df.list.0.5, rep(1500, 20), rep(5000, 20), rep(0.5, 20)), fit_model_M, TRUE)
saveRDS(res.0.5.1M, "res.sim.1M.1500.region.5000.geom.ibdclust2cM.MAF.0.5.boot.keep.ends.rds")
