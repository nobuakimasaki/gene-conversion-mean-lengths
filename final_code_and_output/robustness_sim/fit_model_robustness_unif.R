library(data.table)
library(purrr)
library(parallel)
source("fns_sim.R")
set.seed(26)

# Detect the number of available cores
num_cores <- detectCores()

read.MAF <- function(seed) {
  str = paste0("/projects/browning/brwnlab/sharon/for_nobu/gc_length/sim5_data/sim5_seed", seed,
               "_10Mb_n125000.gtstats")
  df <- read.table(str)
  return(df)
}

# read.MAF <- function(seed) {
#   str = paste0("sim5_seed", seed,
#                "_10Mb_n125000.gtstats")
#   df <- read.table(str)
#   return(df)
# }

boot_MLE_M <- function(x, l_lst, psi_lst, M) {
  ind <- sample(1:length(l_lst), replace = TRUE)
  l_lst_boot <- l_lst[ind]
  psi_lst_boot <- psi_lst[ind]
  
  optim.out.geom <- optim(0.005, neg_log_lik, psi_lst = psi_lst_boot, pL = pL_geom_2M,
                          l_lst = l_lst_boot, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  optim.out.geom2 <- optim(0.005, neg_log_lik, psi_lst = psi_lst_boot, pL = pL_geom2_2M,
                           l_lst = l_lst_boot, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  
  print("finished 1 bootstrap iteration")
  return(c(1/optim.out.geom$par, 2/optim.out.geom2$par))
}

fit_model_M <- function(df, MAF.df, M, region, MAF.ceil) {
  
  # to make indexing easier
  MAF.df$pos <- MAF.df$V2 + 1
  psi <- numeric(max(MAF.df$pos))
  psi[c(MAF.df$pos)] <- 2*MAF.df$V11*(1-MAF.df$V11)
  
  exc <- MAF.df[MAF.df$V11 < 0.05 | MAF.df$V11 > MAF.ceil ,]
  psi[exc$pos] <- 0
  
  # get obs tracts
  df <- dplyr::select(df, -iter)
  df$V1 <- df$V1 + 1
  df$V2 <- df$V2 + 1
  
  tracts_lst <- split(df, seq(nrow(df)))
  psi_lst <- lapply(tracts_lst, est_psi, psi = psi, region = region, length_chrom = length(psi), debias = FALSE)
  l_lst <- lapply(tracts_lst, function(x) {x[2]-x[1]+1})
  
  keep <- which( df$V2 - df$V1 + 1 <= M & df$V2 - df$V1 + 1 > 1)
  
  psi_lst_filt <- psi_lst[keep]
  l_lst_filt <- l_lst[keep]
  
  print(psi_lst_filt %>% unlist())
  print(l_lst_filt %>% unlist())
  
  optim.out.geom <- optim(0.005, neg_log_lik, psi_lst = psi_lst_filt, pL = pL_geom_2M,
                          l_lst = l_lst_filt, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  optim.out.geom2 <- optim(0.005, neg_log_lik, psi_lst = psi_lst_filt, pL = pL_geom2_2M,
                           l_lst = l_lst_filt, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  
  res.boot <- mclapply(1:500, boot_MLE_M, l_lst_filt, psi_lst_filt, M, mc.cores = num_cores)
  
  print("finished 1 replicate")
  print(1/optim.out.geom$par)
  return(list(optim.out.geom, optim.out.geom2, l_lst, psi_lst, res.boot))
}

MAF.chrom.1 <- read.MAF(1)

tracts_unif <- read.csv("sim_tracts_vcf_unif_multiple_iterations.csv")
colnames(tracts_unif) <- c("V1", "V2", "length", "iter")
#tracts_unif <- dplyr::filter(tracts_unif, iter <= 19)
tracts_unif <- dplyr::select(tracts_unif, -length)
tracts_unif_df_list <- split(tracts_unif, tracts_unif$iter)
res_unif <- lapply(tracts_unif_df_list, fit_model_M, MAF.chrom.1, 1500, 5000, 0.5)

saveRDS(res_unif, "res.sim.2M.1500.region.5000.unif.MAF.0.5.boot.keep.ends.rds")

