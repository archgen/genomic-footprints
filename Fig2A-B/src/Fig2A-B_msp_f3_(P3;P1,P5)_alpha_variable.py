#!/usr/bin/env python3
# F3 Statistic Power Simulations
# Figure 2 A-B (varying alpha)
# This Python script simulates admixture events under a simple demographic model using msprime and stdpopsim, computes branch-based 
# f2 and f3 statistics across a range # of admixture proportions (α), saves the results to disk, and generates a bar‐and‐error‐bar 
# plot of the statistics as a function of α.

# Williams & Huber 2025 – Genomic Footprints 

import msprime                    # Coalescent simulator
import tskit                      # Tree sequence toolkit
import matplotlib.pyplot as plt   # Plotting library
import numpy as np                # Numerical computations
import demesdraw                  # Demographic model visualization
import os                         # File system operations
import pandas as pd               # Data manipulation
import seaborn as sns             # Statistical plotting
import stdpopsim                  # Standardized population models

# Output directory for saving figures
outdir = "Users/mkw5910/Documents/admixture_aDNA_review/figures/f3_X_A_B_model/"
# Create the directory if it does not exist
os.makedirs(outdir, exist_ok=True)

# Define a helper function to compute coalescence probability
def probability_of_coalescence(Ne, t):
    """
    Calculate the probability that two lineages coalesce before time t.

    Parameters:
      Ne (int): Effective population size.
      t  (int): Time in generations before present.

    Returns:
      float: Probability of coalescence before time t.
    """
    # Exponential approximation for large Ne
    return 1 - np.exp(-t / (2 * Ne))

# Set up genome properties via stdpopsim
engine  = stdpopsim.get_engine("msprime")
species = stdpopsim.get_species("HomSap")                          # Homo sapiens
demog   = species.get_demographic_model("AncientEurasia_9K19")     # Pre-defined model
contig21 = species.get_contig("chr21", mutation_rate=demog.mutation_rate)
# Total sequence length from recombination map
L  = contig21.recombination_map.position[-1]
# Recombination rate
ro = contig21.recombination_map.rate[0]

# Initialize lists to collect summary statistics
differences             = []
f3_values_mean          = []
f3_values_std           = []
f3_f2_values_mean       = []
f3_f2_values_std        = []
f3_equal_ratio          = []
f2_p3_p1_values_mean    = []
f2_p3_p5_values_mean    = []
f2_p1_p5_values_mean    = []
f2_p3_p1_values_std     = []
f2_p3_p5_values_std     = []
f2_p1_p5_values_std     = []

# Parameters for the simulation
split_gen      = 100       # Generation at which the final split occurs
ne             = 10e2      # Effective population size
alpha          = 0         # Admixture proportion
samples        = 20        # Samples per population
seed           = 77        # RNG seed
num_replicates = 1         # Number of replicates per scenario

# Range of admixture times (from 6 up to split_gen, stepping by 5)
time_values = np.arange(6, split_gen, 5)

# Loop over each admixture time
for time in time_values:
    # Build a demographic model in msprime
    demography = msprime.Demography()
    demography.add_population(name="P1",  initial_size=ne)  # Population P1
    demography.add_population(name="P5",  initial_size=ne)  # Population P5
    demography.add_population(name="P13", initial_size=ne)  # Ancestor P13
    demography.add_population(name="P3",  initial_size=ne)  # Population P3
    demography.add_population(name="P513",initial_size=ne)  # Root P513

    # Introduce admixture: P3 → P5 at generation 5
    demography.add_mass_migration(time=5, source="P3", dest="P5", proportion=alpha)
    # Split P13 into P3 and P1 one generation after admixture
    demography.add_population_split(time=time+1, derived=["P3","P1"], ancestral="P13")
    # Final split of P513 into P13 and P5 at split_gen
    demography.add_population_split(time=split_gen, derived=["P13","P5"], ancestral="P513")

    # Visualize the demographic history with demesdraw
    graph      = demography.to_demes()
    demes_names = [d.name for d in graph.demes]
    palette     = sns.color_palette("colorblind", len(demes_names))
    colours     = dict(zip(demes_names, palette))

    ax = demesdraw.tubes(
        graph,
        colours="dimgrey",   # Tube color
        log_time=True,       # Log-scale time axis
        labels="xticks-mid", # Label placement
        fill=True            # Fill tubes
    )

    # Save the model figure
    filename = os.path.join(
        outdir,
        f"f3(X;A,B)_model_admixG_{time}__splitG_{split_gen}.pdf"
    )
    ax.figure.savefig(filename)

    # Prepare containers for replicate statistics
    f3_replicates       = []
    f3_f2_replicates    = []
    f2_p3_p1_replicates = []
    f2_p3_p5_replicates = []
    f2_p1_p5_replicates = []

    # Simulate ancestry under this demographic model
    replicates = msprime.sim_ancestry(
        {0: samples, 1: samples, 2: samples, 3: samples, 4: samples},  # sample counts
        sequence_length=L,
        demography=demography,
        recombination_rate=ro,
        num_replicates=num_replicates,
        random_seed=seed,
        model=[msprime.DiscreteTimeWrightFisher(duration=20),
               msprime.StandardCoalescent()]
    )

    # Compute f-statistics per replicate
    for ts in replicates:
        # f3(P3; P1, P5)
        f3_value = ts.f3(
            sample_sets=[
                ts.samples(population=3),
                ts.samples(population=1),
                ts.samples(population=0)
            ],
            mode="branch"
        )
        # Pairwise f2 statistics
        f2_p3_p1 = ts.f2(
            sample_sets=[ts.samples(population=3), ts.samples(population=0)],
            mode='branch'
        )
        f2_p3_p5 = ts.f2(
            sample_sets=[ts.samples(population=3), ts.samples(population=1)],
            mode='branch'
        )
        f2_p1_p5 = ts.f2(
            sample_sets=[ts.samples(population=0), ts.samples(population=1)],
            mode='branch'
        )
        # Alternative f3 via f2 values
        f3_f2vals = 0.5 * (f2_p3_p1 + f2_p3_p5 - f2_p1_p5)

        # Store replicate results
        f3_replicates.append(f3_value)
        f3_f2_replicates.append(f3_f2vals)
        f2_p3_p1_replicates.append(f2_p3_p1)
        f2_p3_p5_replicates.append(f2_p3_p5)
        f2_p1_p5_replicates.append(f2_p1_p5)

    # Aggregate replicate means and standard deviations
    f3_values_mean.append(np.mean(f3_replicates))
    f3_values_std.append(np.std(f3_replicates))
    f3_f2_values_mean.append(np.mean(f3_f2_replicates))
    f3_f2_values_std.append(np.std(f3_f2_replicates))
    f2_p3_p1_values_mean.append(np.mean(f2_p3_p1_replicates))
    f2_p3_p5_values_mean.append(np.mean(f2_p3_p5_replicates))
    f2_p1_p5_values_mean.append(np.mean(f2_p1_p5_replicates))
    f2_p3_p1_values_std.append(np.std(f2_p3_p1_replicates))
    f2_p3_p5_values_std.append(np.std(f2_p3_p5_replicates))
    f2_p1_p5_values_std.append(np.std(f2_p1_p5_replicates))

    # Calculate theoretical f3_equal ratio using coalescence probability
    tr   = split_gen
    t1   = time
    prob = probability_of_coalescence(ne, t1)
    f3_equal = ((1 / (1 - prob)) * (t1 / tr)) - 2 * alpha * (1 - alpha)
    f3_equal_ratio.append(f3_equal)

    # Track difference between split and admixture times
    differences.append(split_gen - time)

## PLOTTING

import matplotlib.pyplot as plt  # Redundant import for plotting
import numpy as np               # Redundant import for numeric ops

# Assemble results into a DataFrame
data = {
    'f2_p3_p1':     f2_p3_p1_values_mean,
    'f2_p3_p1_sd':  f2_p3_p1_values_std,
    'f2_p3_p5':     f2_p3_p5_values_mean,
    'f2_p3_p5_sd':  f2_p3_p5_values_std,
    'f2_p1_p5':     f2_p1_p5_values_mean,
    'f2_p1_p5_sd':  f2_p1_p5_values_std,
    'f3(pX; pA, pB)': 0.5 * (
          np.array(f2_p3_p1_values_mean) +
          np.array(f2_p3_p5_values_mean) -
          np.array(f2_p1_p5_values_mean)
    ),
    'f3_se':        f3_values_std,
    'Diff_splitgen_admixgen': time_values - 5
}
df = pd.DataFrame(data)
# Preview the compiled table
print(df)

# Create combined bar + error-bar line plot
plt.clf()
sns.set(style="whitegrid", font_scale=1.3)
cb = sns.color_palette("colorblind", 4)
col_f3, col_f2_p3p1, col_f2_p3p5, col_f2_p1p5 = cb

fig, ax = plt.subplots(figsize=(12, 6))

# Bar for f3 with error bars
ax.bar(
    df['Diff_splitgen_admixgen'],
    df['f3(pX; pA, pB)'],
    width=0.8,
    color=col_f3,
    yerr=df['f3_se'],
    capsize=5,
    label='f3(P3; P1, P5) as 1/2( f2(P3,P1) + f2(P3,P5) - f2(P1,P5) )'
)

# Line plots with error bars for each f2
ax.errorbar(
    df['Diff_splitgen_admixgen'],
    df['f2_p3_p1'],
    yerr=df['f2_p3_p1_sd'],
    color=col_f2_p3p1,
    marker='o',
    linewidth=2,
    fmt='-o',
    capsize=5,
    label='f2(P3,P1)'
)
ax.errorbar(
    df['Diff_splitgen_admixgen'],
    df['f2_p3_p5'],
    yerr=df['f2_p3_p5_sd'],
    color=col_f2_p3p5,
    marker='o',
    linewidth=2,
    fmt='-o',
    capsize=5,
    label='f2(P3,P5)'
)
ax.errorbar(
    df['Diff_splitgen_admixgen'],
    df['f2_p1_p5'],
    yerr=df['f2_p1_p5_sd'],
    color=col_f2_p1p5,
    marker='o',
    linewidth=2,
    fmt='-o',
    capsize=5,
    label='f2(P1,P5)'
)

# Horizontal zero line
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# Label axes
ax.set_xlabel(
    f'Difference in (P1, P3) split time and admixture generation ((α = {alpha}) P5 → P3)',
    fontsize=18
)
ax.set_ylabel('f-statistic estimate', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=16)

# Adjust layout and y-axis limits
plt.tight_layout()
ax.set_ylim(-50, 250)

# Extract legend to a separate figure
fig_legend, ax_legend = plt.subplots(figsize=(8, 1.5))
ax_legend.axis('off')
handles, labels = ax.get_legend_handles_labels()
ax_legend.legend(handles, labels, fontsize=12, loc='center', frameon=False)

plt.tight_layout()
plt.show()
