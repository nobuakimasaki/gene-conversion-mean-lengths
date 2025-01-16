This directory contains the code and data used to fit the model on the UK Biobank sequence data. The multi-individual IBD method is used to first detect observed tracts from the UK Biobank sequence data. 

The "data" directory contains the detected tracts using the multi-individual IBD method (by chromosome).

In the "hotspot_detection" directory, we have the deCode 2019 genetic maps for each chromosome. In this directory, "detection.R" is first used to calculate recombination rates between pairs of markers in the genetic maps. This is saved in "res_recombination_rates.rds". Then, "detection_and_analysis.Rmd" is used to detect recombination hotspots, which is saved in "hotspots.RDS". 

"fit_model_M_UK_stratified.R" is then used to fit the model on the detected tracts. We fit the model on all detected tracts, and also stratifying by the hotspots in "hotspots.RDS". "fns_sim.R" contains functions used to fit the model. "fit_model_M_UK_stratified.R" also uses minor allele frequencies from all 22 autosomes, but we cannot include this data. Results are then saved in "res.UK_Biobank.1M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.grouped.keep.ends.hotspot.rds", "res.UK_Biobank.1M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.grouped.keep.ends.not.hotspot.rds", and "res.UK_Biobank.1M.1500.region.5000.ibdclust2cM.MAF.0.5.boot.grouped.keep.ends.rds".

These results are analyzed in "hotspot_detection/detection_and_analysis.Rmd". 


