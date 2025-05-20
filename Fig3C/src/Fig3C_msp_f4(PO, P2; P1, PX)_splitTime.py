#!/usr/bin/env python3
# F4 Statistic Power Simulations
# Fig 3C
# This script simulates admixture events (P3 ← α·P1 + (1–α)·P5) under a simple demographic model 
# using msprime and stdpopsim, computes branch‐based f2 and f3 statistics across a range of admixture proportions (α), 
# saves the results, and visualizes how these statistics vary with α.

# Williams & Huber 2025 – Genomic Footprints 

import os                                    # file system operations
import numpy as np                           # numerical computing
import pandas as pd                          # data manipulation
import matplotlib.pyplot as plt              # plotting
import seaborn as sns                        # statistical plotting
import msprime                               # coalescent simulations
import tskit                                 # tree sequence utilities
import demesdraw                             # demographic model plotting
import stdpopsim                             # standardized population models

# Set output directory (adjust path as needed) and create it if it doesn't exist
outdir = "/Users/mkw5910/Documents/admixture_aDNA_review/figures/f3_X_A_B_model/"
os.makedirs(outdir, exist_ok=True)

# Function: probability of coalescence
def probability_of_coalescence(Ne, t):
    """
    Calculate the probability that two lineages coalesce before time t.
    
    Parameters:
      Ne (int): Effective population size (assumed constant).
      t (int): Time before present (in generations).
      
    Returns:
      float: The probability of coalescence before time t.
    """
    return 1 - np.exp(-t / (2 * Ne))

# Load genome properties using stdpopsim
engine = stdpopsim.get_engine("msprime")                 # simulation engine
species = stdpopsim.get_species("HomSap")                # species: Homo sapiens
demog = species.get_demographic_model("AncientEurasia_9K19")
contig21 = species.get_contig("chr21", mutation_rate=demog.mutation_rate)
L = contig21.recombination_map.position[-1]              # sequence length
ro = contig21.recombination_map.rate[0]                  # recombination rate

# Initialize lists to store summary statistics across split times
f4_power_iter       = []
differences         = []
f4_values_mean      = []
f4_values_std       = []
f4_f2_values_mean   = []
f4_f2_values_std    = []

f2_PO_P1_values_mean = []
f2_P2_PX_values_mean = []
f2_PO_PX_values_mean = []
f2_P1_P2_values_mean = []

f2_PO_P1_values_std  = []
f2_P2_PX_values_std  = []
f2_PO_PX_values_std  = []
f2_P1_P2_values_std  = []

# Simulation parameters
split_gen      = 100        # time of P1P2 vs PO split
ne             = 1000       # effective population size
alpha          = 0.25       # admixture proportion
samples        = 20         # number of samples per population
seed           = 77         # random seed for reproducibility
num_replicates = 1          # number of simulation replicates
admix_gen      = 5          # admixture generation

# Define range of split times to iterate over
time_values = np.arange(6, split_gen, 5)

# Loop over each split time
for time in time_values:
    # Build the demographic model
    demography = msprime.Demography()
    demography.add_population(name="P1", initial_size=ne)      # pop index 0
    demography.add_population(name="P2", initial_size=ne)      # pop index 1
    demography.add_population(name="P1P2", initial_size=ne)    # pop index 2
    demography.add_population(name="PX", initial_size=ne)      # pop index 3
    demography.add_population(name="PO", initial_size=ne)      # pop index 4
    demography.add_population(name="PANC", initial_size=ne)    # pop index 5

    # Add admixture event: P1 + P2 → PX at admix_gen
    demography.add_admixture(
        time=admix_gen,
        derived="PX",
        ancestral=["P1", "P2"],
        proportions=[alpha, (1 - alpha)]
    )
    # Split P1 and P2 from P1P2 at 'time'
    demography.add_population_split(
        time=time,
        derived=["P1", "P2"],
        ancestral="P1P2"
    )
    # Split P1P2 and PO from PANC at split_gen
    demography.add_population_split(
        time=split_gen,
        derived=["P1P2", "PO"],
        ancestral="PANC"
    )

    # Visualize and save the demographic model
    graph = demography.to_demes()
    ax_model = demesdraw.tubes(
        graph,
        colours="dimgrey",
        log_time=True,
        labels="xticks-mid",
        fill=True
    )
    model_filename = os.path.join(
        outdir,
        f"f4(PO,P2;P1,PX)_model_admixG_{time}_splitG_{split_gen}.pdf"
    )
    ax_model.figure.savefig(model_filename)    # save PDF
    plt.close(ax_model.figure)                # free memory

    # Initialize per-replicate lists
    f4_replicates       = []
    f4_f2_replicates    = []
    f2_PO_P1_replicates = []
    f2_P2_PX_replicates = []
    f2_PO_PX_replicates = []
    f2_P1_P2_replicates = []

    # Run coalescent simulations
    replicates = msprime.sim_ancestry(
        samples={0: samples, 1: samples, 2: samples, 3: samples, 4: samples},
        sequence_length=L,
        demography=demography,
        recombination_rate=ro,
        num_replicates=num_replicates,
        random_seed=seed,
        model=[msprime.DiscreteTimeWrightFisher(duration=20),
               msprime.StandardCoalescent()]
    )
    
    # Process each tree sequence replicate
    for ts in replicates:
        # Compute f4 statistic directly
        f4_val = ts.f4(
            sample_sets=[
                ts.samples(population=4),  # PO
                ts.samples(population=1),  # P2
                ts.samples(population=0),  # P1
                ts.samples(population=3)   # PX
            ],
            mode="branch"
        )
        # Compute pairwise f2 statistics
        f2_PO_P1 = ts.f2(
            sample_sets=[ts.samples(population=4), ts.samples(population=0)],
            mode="branch"
        )
        f2_P2_PX = ts.f2(
            sample_sets=[ts.samples(population=1), ts.samples(population=3)],
            mode="branch"
        )
        f2_PO_PX = ts.f2(
            sample_sets=[ts.samples(population=4), ts.samples(population=3)],
            mode="branch"
        )
        f2_P1_P2 = ts.f2(
            sample_sets=[ts.samples(population=0), ts.samples(population=1)],
            mode="branch"
        )
        
        # Derive f4 from f2 values
        f4_from_f2 = 0.5 * (
            f2_PO_PX + f2_P1_P2 - f2_PO_P1 - f2_P2_PX
        )
        
        # Append replicate values
        f4_replicates.append(f4_val)
        f4_f2_replicates.append(f4_from_f2)
        f2_PO_P1_replicates.append(f2_PO_P1)
        f2_P2_PX_replicates.append(f2_P2_PX)
        f2_PO_PX_replicates.append(f2_PO_PX)
        f2_P1_P2_replicates.append(f2_P1_P2)
    
    # Compute mean and standard deviation for this split time
    f4_values_mean.append(np.mean(f4_replicates))
    f4_values_std.append(np.std(f4_replicates))
    f4_f2_values_mean.append(np.mean(f4_f2_replicates))
    f4_f2_values_std.append(np.std(f4_f2_replicates))
    
    f2_PO_P1_values_mean.append(np.mean(f2_PO_P1_replicates))
    f2_P2_PX_values_mean.append(np.mean(f2_P2_PX_replicates))
    f2_PO_PX_values_mean.append(np.mean(f2_PO_PX_replicates))
    f2_P1_P2_values_mean.append(np.mean(f2_P1_P2_replicates))
    
    f2_PO_P1_values_std.append(np.std(f2_PO_P1_replicates))
    f2_P2_PX_values_std.append(np.std(f2_P2_PX_replicates))
    f2_PO_PX_values_std.append(np.std(f2_PO_PX_replicates))
    f2_P1_P2_values_std.append(np.std(f2_P1_P2_replicates))
    
    # Track difference between split and admixture times
    differences.append(split_gen - time)
    # f4 power parameter: (1 – α) * (time – admix_gen)
    f4_power_iter.append((1 - alpha) * (time - admix_gen))

# Build summary DataFrame for plotting
data = {
    'f2_PO_P1':       f2_PO_P1_values_mean,
    'f2_PO_P1_sd':    f2_PO_P1_values_std,
    'f2_P2_PX':       f2_P2_PX_values_mean,
    'f2_P2_PX_sd':    f2_P2_PX_values_std,
    'f2_PO_PX':       f2_PO_PX_values_mean,
    'f2_PO_PX_sd':    f2_PO_PX_values_std,
    'f2_P1_P2':       f2_P1_P2_values_mean,
    'f2_P1_P2_sd':    f2_P1_P2_values_std,
    'f4(PO,P2;P1,PX)': 0.5 * (
        np.array(f2_PO_PX_values_mean) + 
        np.array(f2_P1_P2_values_mean) -
        np.array(f2_P2_PX_values_mean) -
        np.array(f2_PO_P1_values_mean)
    ),
    'f4_se':          f4_values_std,
    'Diff_admix_time': abs(admix_gen - time_values),
    'f4_power_iter':  f4_power_iter
}
df = pd.DataFrame(data)
print(df)

# ----- Plotting: f4 statistic and f2 components vs. admixture time difference -----

sns.set(style="whitegrid", font_scale=1.3)          # set plotting style
cb  = sns.color_palette("colorblind", 5)            # colorblind palette
col_f4, col_PO_P1, col_PO_PX, col_P2_PX, col_P1_P2 = cb

# Simple bar + errorbar plot
fig, ax = plt.subplots(figsize=(12, 8))
ax.bar(
    df['Diff_admix_time'],
    df['f4(PO,P2;P1,PX)'],
    width=0.8,
    color=col_f4,
    yerr=df['f4_se'],
    capsize=5,
    label='f4(PO,P2;P1,PX)'
)
ax.errorbar(
    df['Diff_admix_time'], df['f2_PO_P1'],
    yerr=df['f2_PO_P1_sd'], color=col_PO_P1,
    marker='o', linewidth=2, fmt='-o', capsize=5,
    label='f2(PO,P1)'
)
ax.errorbar(
    df['Diff_admix_time'], df['f2_P2_PX'],
    yerr=df['f2_P2_PX_sd'], color=col_P2_PX,
    marker='o', linewidth=2, fmt='-o', capsize=5,
    label='f2(P2,PX)'
)
ax.errorbar(
    df['Diff_admix_time'], df['f2_PO_PX'],
    yerr=df['f2_PO_PX_sd'], color=col_PO_PX,
    marker='o', linewidth=2, fmt='-o', capsize=5,
    label='f2(PO,PX)'
)
ax.errorbar(
    df['Diff_admix_time'], df['f2_P1_P2'],
    yerr=df['f2_P1_P2_sd'], color=col_P1_P2,
    marker='o', linewidth=2, fmt='-o', capsize=5,
    label='f2(P1,P2)'
)
ax.scatter(
    df['Diff_admix_time'], df['f4_power_iter'],
    color='black', marker='*', s=100,
    label='(1–α)*(time – admix_gen)'
)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# Axis labels and formatting
ax.set_xlabel(
    f'Difference in split vs. admixture generation (α={alpha})',
    fontsize=19
)
ax.set_ylabel('f-statistic estimate', fontsize=19)
ax.set_ylim(-20, 250)
ax.tick_params(axis='both', which='major', labelsize=19)
ax.grid(True)

plt.tight_layout()
plt.show()

# BIG Legend Plot
fig, (ax, ax_legend) = plt.subplots(
    nrows=2,
    figsize=(12, 20),
    gridspec_kw={'height_ratios': [2, 3]}
)
# (repeat plotting commands on ax...)
# [omitted here for brevity; same as above but with inverted x-axis and custom legend]
# plt.tight_layout() and plt.show() at end
