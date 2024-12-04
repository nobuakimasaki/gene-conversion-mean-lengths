library(data.table)
library(purrr)
library(parallel)
library(dplyr)
source("fns_sim.R")
#chr <- commandArgs(trailingOnly = TRUE)
set.seed(26)

check_tract_in_hotspot <- function(tract, hotspots) {
  tract_start <- tract[[1]]
  tract_end <- tract[[2]]
  
  hotspots$max_start <- pmax(hotspots$first_pos, tract_start)
  hotspots$min_end <- pmin(hotspots$last_pos, tract_end)
  hotspots$overlap <- hotspots$max_start <= hotspots$min_end
  
  # Return TRUE if dataset not empty, otherwise FALSE
  return(sum(hotspots$overlap) > 0)
}

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

tracts.df.chr.21 <- tracts.df.0.5.lst[[21]]
# Remove singleton tracts and tracts greater than ceiling
keep <- which( tracts.df.chr.21$V2 - tracts.df.chr.21$V1 + 1 <= 1500 & tracts.df.chr.21$V2 - tracts.df.chr.21$V1 + 1 > 1)
tracts.df.chr.21 <- tracts.df.chr.21[keep, ]

# Obtain full tract and psi list
tracts_lst.chr.21 <- split(tracts.df.chr.21, seq(nrow(tracts.df.chr.21)))
hotspots.chr.21 <- hotspots.chrs[[21]] %>% filter(hotspot == TRUE)

inside_hotspot <- lapply(tracts_lst.chr.21, check_tract_in_hotspot, hotspots.chr.21)

saveRDS(list(tracts_lst.chr.21, inside_hotspot), "inside_hotspot.chr.21.rds")