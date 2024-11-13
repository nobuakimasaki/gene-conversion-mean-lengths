import allel
import time
import numpy as np
import pandas as pd
import csv
import multiprocessing
from functools import partial
import random

random.seed(27)

print("n cores:")
print(multiprocessing.cpu_count())

def sim_tracts(N_gene_conv, min_pos, max_pos, dist_N, mean = 300):
    # Generate start positions
    start = np.random.choice(range(min_pos, max_pos-2000), N_gene_conv, replace=True)

    # Generate lengths based on the specified distribution
    if dist_N == "geom":
        lengths = np.random.geometric(1/mean, N_gene_conv)
    elif dist_N == "geom2":
        lengths = np.random.geometric(2/mean, N_gene_conv) + np.random.geometric(2/mean, N_gene_conv)
    elif dist_N == "unif":
        lengths = np.random.uniform(low = 1, high = mean*2-1, size = N_gene_conv)
    elif dist_N == "negbinom":
        lengths = np.random.negative_binomial(3, 3/(3+mean), size = N_gene_conv)
    else:
        raise ValueError("Invalid distribution type.")

    # Calculate end positions
    end = start + lengths - 1
    # Combine start and end positions into a single array
    start_end = np.column_stack((start, end))
    # Filter out rows where the end position is greater than max_pos
    start_end = start_end[start_end[:, 1] < max_pos]
    # Split the start_end array into a list of arrays (one per tract)
    tracts = [start_end[i] for i in range(start_end.shape[0])]

    return tracts

def sim_gene_conv(actual, idx, positions, genotypes):
    # print("inside sim_gene_conv function")
    # print("actual: ")
    # print(actual)
    # print("positions: ")
    # print(positions)

    geno_series = genotypes[:, idx]
    # Convert geno_series to a NumPy array for easier manipulation
    geno_array = np.array(geno_series)
    # Create a new list with 1 where geno_series element is [0, 0] and 0 otherwise
    heterozygous = [1 if np.array_equal(geno, [1, 0]) or np.array_equal(geno, [0, 1]) else 0 for geno in geno_array]
    heterozygous = np.array(heterozygous)

    # print("heterozygous: ")
    # print(heterozygous)
    # # Define the start and end of the actual tract
    start, end = actual
    # print("start: ")
    # print(start)
    # print("end: ")
    # print(end)
    # Identify the indices of the positions within the tract
    tract_ind = np.where((positions >= start) & (positions <= end))[0]
    tract_ind = tract_ind.astype(int)
    print("tract_ind: ")
    print(tract_ind)

    if len(tract_ind) == 0:
        # print("L: ")
        # print(0)
        return [0, 0, 0]

    # Extract the heterozygous markers within the tract
    tract = heterozygous[tract_ind]
    # print("tract: ")
    # print(tract)
    # print("positions[tract_ind]: ")
    # print(positions[tract_ind])

    # Calculate the length of the gene conversion tract
    if np.sum(tract) == 0:
        result = [0, 0, 0]
    elif np.sum(tract) == 1:
        obs_start = positions[tract_ind[np.min(np.where(tract == 1))]]
        obs_end = positions[tract_ind[np.max(np.where(tract == 1))]]
        L = obs_end - obs_start + 1
        result = [obs_start, obs_end, L]
    else:
        # Get the index of the leftmost and rightmost markers and locate the positions
        obs_start = positions[tract_ind[np.min(np.where(tract == 1))]]
        obs_end = positions[tract_ind[np.max(np.where(tract == 1))]]
        L = obs_end - obs_start + 1
        result = [obs_start, obs_end, L]

    print("result: ")
    print(result)
    return result

# Define the file path
file_path = "/projects/browning/brwnlab/sharon/for_nobu/gc_length/sim5_data/sim5_seed1_10Mb_n125000.gtstats"
# Read the table into a pandas DataFrame
maf_df = pd.read_table(file_path, header=None)  # Assuming the file has no header
# Filter the DataFrame and select column V2 where V11 < 0.05
keep = maf_df[maf_df.iloc[:, 10] >= 0.05].iloc[:, 1].astype(int).tolist()
print("keep: ")
print(keep)
###### So far, we've defined the functions and the MAF file that will be used to generate the tracts.
###### We next load in the genotypes for individuals and actually simulate the tracts.

# Path to your compressed VCF file
vcf_file = "/projects/browning/brwnlab/sharon/for_nobu/gc_length/sim5_vcfs/sim5_seed1_10Mb_n125000_err0.0002phased_del1.vcf.gz"
# vcf_file = "example.vcf"
# Open the compressed VCF file for reading
callset = allel.read_vcf(vcf_file)
print("callset: ")
print(sorted(callset.keys()))
# Extract the samples, genotype calls, and variant positions
samples = callset['samples']
genotype_calls = callset['calldata/GT']
variant_positions = callset['variants/POS']
variant_positions = variant_positions.astype(int)

filtered_indices = [pos in keep for pos in variant_positions]
filtered_positions = variant_positions[filtered_indices]
filtered_genotypes = genotype_calls[filtered_indices]

print(filtered_genotypes[:5, :5])

# Specify the number of individuals you want to sample
N = 10**5
num_iterations = 100  # Specify the number of iterations

all_data = []  # List to store the concatenated results

for iteration in range(num_iterations):
    print(iteration)
    
    # Simulate gene conversion tracts
    gene_conversion_tracts = sim_tracts(N, min(filtered_positions), max(filtered_positions), "negbinom")
    
    # Define function with prespecified positions
    fixed_positions_sim_gene_conv = partial(sim_gene_conv, positions=filtered_positions, genotypes=filtered_genotypes)
    
    data = []
    for i in range(N):
        # Randomly sample one individual and obtain genotype data
        idx = np.random.choice(filtered_genotypes.shape[1])
        # Get one gene conversion tract
        actual = gene_conversion_tracts[i]
        data.append((actual, idx))

    # Use multiprocessing pool to run the simulation
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1)
    L = pool.starmap(fixed_positions_sim_gene_conv, data)
    pool.close()
    pool.join()

    # Add the iteration number as the fourth element in each sublist
    L_with_iteration = [lst + [iteration] for lst in L if lst[-1] != 0]  # Filter out sublists where the last element is 0
    
    # Concatenate the results from this iteration
    all_data.extend(L_with_iteration)

# Write the concatenated results to a CSV file
with open('sim_tracts_vcf_negbinom_multiple_iterations.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(all_data)
