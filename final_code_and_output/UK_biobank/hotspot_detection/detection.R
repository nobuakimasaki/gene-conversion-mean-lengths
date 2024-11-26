library(dplyr)

calc_cM_per_bp <- function(df) {
  first_row <- 1
  n <- nrow(df)
  results <- list()
  j <- 0
  # Loop through each row
  while (j < n & first_row < n) {
    # Loop to find the next row that is at least 2000 (2 kb) apart
    for (j in (first_row + 1):n) {
      # Check if the distance is at least 2000
      if (df$V4[j] - df$V4[first_row] >= 2000) {
        # Calculate the difference in V3 divided by the difference in V4
        result <- (df$V3[j] - df$V3[first_row]) / (df$V4[j]/10^6 - df$V4[first_row]/10^6)
        # Store the result
        results[[length(results) + 1]] <- c(result, first_row, j, df$V4[first_row], df$V4[j])
        first_row <- j
        break
      }
    }
  }
  res.df <- results %>% unlist() %>% matrix(byrow = TRUE, ncol = 5) %>% as.data.frame()
  res.df$center <- (res.df$V4 + res.df$V5)/2
  
  background.rate <- (df$V3[nrow(df)] - df$V3[1]) / (df$V4[nrow(df)]/10^6 - df$V4[1]/10^6)
  
  return(list(res.df, background.rate))
}

# read_map_file <- function(chr) {
#   str = paste0("/projects/browning/maps/decode.2019.b38/decode2019.chrchr", chr,
#                ".GRCh38.map")
#   df <- read.table(str)
#   return(df)
# }

read_map_file <- function(chr) {
  str = paste0("decode2019.chrchr", chr,
               ".GRCh38.map")
  df <- read.table(str)
  return(df)
}

map_list <- list()

for (i in 1:22) {
  # Construct the chromosome name as a string
  chr_name <- paste0("chr", i, ".map")
  # Read the map file and store it in the list
  map_list[[chr_name]] <- read_map_file(as.character(i))
}

res <- lapply(map_list, calc_cM_per_bp)
saveRDS(res, "res_recombination_rates.rds")

print(head(res[[1]][[1]]))
#ggplot(res.df, aes(center, V1)) + geom_point() 