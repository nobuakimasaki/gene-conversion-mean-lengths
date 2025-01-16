This directory contains the code and data used to run the simulation in the main paper. In this simulation, gene conversion events are simulated in msprime and gene conversion tracts are subsequently detected using the multi-individual IBD method. 

The "data" directory contains the tract indices found using the multi-individual IBD method, and the minor allele frequencies of each region simulated using msprime.

"fit_model_M_sim_keep_ends.R" uses these tract indices and minor allele frequencies to estimate the mean length of gene conversion tracts for each region. These results are saved in "res.sim.2M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.keep.ends.rds". "fns_sim.R" contains functions loaded and used in "fit_model_M_sim_keep_ends.R".

"res_sim.Rmd" analyzes the results saved in "res.sim.2M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.keep.ends.rds".
