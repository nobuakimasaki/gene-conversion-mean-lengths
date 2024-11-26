library(data.table)
library(purrr)
library(parallel)
library(dplyr)
source("fns_sim.R")
#chr <- commandArgs(trailingOnly = TRUE)
set.seed(26)

# Detect the number of available cores
num_cores <- 20

read.tracts <- function(chr) {
  str <- paste0("/projects/browning/brwnlab/sharon/for_nobu/gc_length/ukbiobank/chr", chr, ".ibdclust2cM_trim1_combinedoffsets_v6.inf_obs_tracts2")
  df <- read.table(str)
  return(df)
}

read.MAF <- function(chr) {
  str <- paste0("/projects/browning/brwnlab/sharon/for_nobu/gc_length/ukbiobank/chr", chr, ".allregions.pmaf.gz")
  df <- read.table(str)
  return(df)
}

tracts.df.0.5.lst <- lapply(1:22, read.tracts)
MAF.df.lst <- lapply(1:22, read.MAF)
hotspots.chrs <- readRDS("hotspots.RDS")

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

check_tract_in_hotspot <- function(tract, hotspots) {
  tract_start <- tract[[1]]
  tract_end <- tract[[2]]
  
  hotspots$max_start <- pmax(hotspots$first_pos, tract_start)
  hotspots$min_end <- pmin(hotspots$last_pos, tract_end)
  hotspots$overlap <- hotspots$max_start <= hotspots$min_end
  
  # Return TRUE if dataset not empty, otherwise FALSE
  return(sum(hotspots$overlap) > 0)
}

### For the UK Biobank data, we want to group all of the tract lengths and psi values from the 22 chromosomes. So we define a function for this.
get_l_psi <- function(MAF.df, df, hotspots, M, region, MAF.ceil) {
  
  psi <- numeric(max(MAF.df$V2))
  psi[c(MAF.df$V2)] <- 2*MAF.df$V3*(1-MAF.df$V3)
  
  # Set psi value to 0 when MAF is less than 0.05 or greater than MAF ceiling
  exc <- MAF.df[MAF.df$V3 < 0.05 | MAF.df$V3 > MAF.ceil ,]
  psi[exc$V2] <- 0
  
  # Remove singleton tracts and tracts greater than ceiling
  keep <- which( df$V2 - df$V1 + 1 <= M & df$V2 - df$V1 + 1 > 1)
  df <- df[keep, ]
  
  print(nrow(df))
  
  # Obtain full tract and psi list
  tracts_lst <- split(df, seq(nrow(df)))
  psi_lst <- lapply(tracts_lst, est_psi, psi = psi, region = region, length_chrom = length(psi), debias = FALSE)
  
  # Get recombination hotspots
  hotspots <- hotspots %>% filter(hotspot == TRUE)
  print("hotspots")
  print(hotspots)
  
  # Get a vector of booleans representing whether tract is in hotspot
  inside_hotspot <- lapply(tracts_lst, check_tract_in_hotspot, hotspots = hotspots) %>% unlist()
  
  # Filter tracts by whether they are in hotspot
  tracts_lst_hotspot <- tracts_lst[inside_hotspot]
  tracts_lst_not_hotspot <- tracts_lst[!inside_hotspot]
  
  # Filter psi by whether corresponding tract is in hotspot
  psi_lst_hotspot <- psi_lst[inside_hotspot]
  psi_lst_not_hotspot <- psi_lst[!inside_hotspot]
  
  l_lst <- lapply(tracts_lst, function(x) {x[[2]]-x[[1]]+1}) 
  l_lst_hotspot <- lapply(tracts_lst_hotspot, function(x) {x[[2]]-x[[1]]+1}) 
  l_lst_not_hotspot <- lapply(tracts_lst_not_hotspot, function(x) {x[[2]]-x[[1]]+1})
  
  print("finished processing data for 1 chromosome")
  
  return(list(l_lst_hotspot, l_lst_not_hotspot, l_lst, psi_lst_hotspot, psi_lst_not_hotspot, psi_lst))
}

l_psi_lst <- pmap(list(MAF.df.lst, tracts.df.0.5.lst, hotspots.chrs), get_l_psi, M = 1500, region = 5000, MAF.ceil = 0.5)

l_lst_hotspot <- lapply(l_psi_lst, function(x) {x[[1]]}) %>% unlist()
l_lst_not_hotspot <- lapply(l_psi_lst, function(x) {x[[2]]}) %>% unlist()
l_lst <- lapply(l_psi_lst, function(x) {x[[3]]}) %>% unlist()

psi_lst_hotspot <- lapply(l_psi_lst, function(x) {x[[4]]}) %>% unlist()
psi_lst_not_hotspot <- lapply(l_psi_lst, function(x) {x[[5]]}) %>% unlist()
psi_lst <- lapply(l_psi_lst, function(x) {x[[6]]}) %>% unlist()

print("l_lst_hotspot")
print(l_lst_hotspot)

print("l_lst_not_hotspot")
print(l_lst_not_hotspot)

print("psi_lst_hotspot")
print(psi_lst_hotspot)

print("psi_lst_not_hotspot")
print(psi_lst_not_hotspot)

fit_model_M <- function(l_lst, psi_lst, M) {
  
  ### maximum likelihood
  optim.out.geom <- optim(0.005, neg_log_lik, psi_lst = psi_lst, pL = pL_geom_2M,
                          l_lst = l_lst, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  optim.out.geom2 <- optim(0.005, neg_log_lik, psi_lst = psi_lst, pL = pL_geom2_2M,
                           l_lst = l_lst, M = M, lower = 0.0001, upper = 0.05, method = "Brent")
  
  res.boot <- mclapply(1:500, boot_MLE_M, l_lst, psi_lst, M, mc.cores = num_cores)
  
  print("finished 1 seed")
  print(1/optim.out.geom$par)
  return(list(optim.out.geom, optim.out.geom2, l_lst, psi_lst, res.boot))
}

res.0.5.1M.hotspot <- fit_model_M(l_lst_hotspot, psi_lst_hotspot, 1500)
res.0.5.1M.not.hotspot <- fit_model_M(l_lst_not_hotspot, psi_lst_not_hotspot, 1500)
res.0.5.1M <- fit_model_M(l_lst, psi_lst, 1500)

saveRDS(res.0.5.1M.hotspot, "res.UK_Biobank.1M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.grouped.keep.ends.hotspot.rds")
saveRDS(res.0.5.1M.not.hotspot, "res.UK_Biobank.1M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.grouped.keep.ends.not.hotspot.rds")
saveRDS(res.0.5.1M, "res.UK_Biobank.1M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.grouped.keep.ends.rds")