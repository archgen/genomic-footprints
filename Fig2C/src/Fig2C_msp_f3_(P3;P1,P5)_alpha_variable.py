#!/usr/bin/env python3
# F3 Statistic Power Simulations
# Fig 2C
# This Python script simulates admixture events under a simple demographic model using msprime and stdpopsim, computes branch-based 
# f2 and f3 statistics across a range # of admixture proportions (α), saves the results to disk, and generates a bar‐and‐error‐bar 
# plot of the statistics as a function of α.

# Williams & Huber 2025 – Genomic Footprints 

"""
f3_admixture_simulation.py

Simulate admixture events (P3 ← α·P1 + (1−α)·P5), compute branch-based f2/f3 statistics,
and visualize how these statistics vary with admixture proportion α.
"""

import os
import msprime           # for coalescent simulations
import tskit             # for tree sequence statistics
import stdpopsim         # to load realistic demographic models
import numpy as np       # numerical computing
import pandas as pd      # tabular data handling
import matplotlib.pyplot as plt
import seaborn as sns    # enhanced plotting
import demesdraw         # (optional) drawing demographies


# -----------------------------------------------------------------------------
# 1. Output directory setup
# -----------------------------------------------------------------------------
outdir = "figures/f3_admixture/"           # relative path to save figures
os.makedirs(outdir, exist_ok=True)         # create dir if needed

# -----------------------------------------------------------------------------
# 2. Utility function
# -----------------------------------------------------------------------------
def probability_of_coalescence(Ne, t):
    """
    Probability that two lineages coalesce before time t generations ago,
    under a constant Ne model.
    """
    return 1 - np.exp(-t / (2 * Ne))


# -----------------------------------------------------------------------------
# 3. Load genome properties & demographic model via stdpopsim
# -----------------------------------------------------------------------------
species = stdpopsim.get_species("HomSap")              # Homo sapiens
demog   = species.get_demographic_model("AncientEurasia_9K19")
contig  = species.get_contig("chr21", mutation_rate=demog.mutation_rate)
L       = contig.recombination_map.position[-1]        # total sequence length
ro      = contig.recombination_map.rate[0]             # uniform recomb. rate

# -----------------------------------------------------------------------------
# 4. Simulation parameters
# -----------------------------------------------------------------------------
split_gen      = 100       # final split time of P13 & P5 from ancestor
XA_split_gen   = split_gen - 90  # time of P3/P1 split
admix_gen      = 5         # generation of admixture event
ne             = int(1e3)  # effective population size
samples        = 20        # samples per population
seed           = 77        # RNG seed
num_reps       = 1         # replicates per admixture proportion

# range of admixture proportions α ∈ [0.0, 0.5] in steps of 0.01
admix_values = np.arange(0, 0.51, 0.01)

# -----------------------------------------------------------------------------
# 5. Prepare storage for summary statistics
# -----------------------------------------------------------------------------
results = {
    "alpha":           [],
    "f2_p3_p1":        [], "sd_f2_p3_p1":        [],
    "f2_p3_p5":        [], "sd_f2_p3_p5":        [],
    "f2_p1_p5":        [], "sd_f2_p1_p5":        [],
    "f3_branch":       [], "sd_f3_branch":       [],
}

# -----------------------------------------------------------------------------
# 6. Main simulation loop
# -----------------------------------------------------------------------------
for alpha in admix_values:
    # 6.1 Build msprime demography
    demography = msprime.Demography()
    demography.add_population(name="P1",   initial_size=ne)
    demography.add_population(name="P5",   initial_size=ne)
    demography.add_population(name="P13",  initial_size=ne)
    demography.add_population(name="P3",   initial_size=ne)
    demography.add_population(name="P513", initial_size=ne)

    # Admixture at time=admix_gen: P3 ← α·P1 + (1−α)·P5
    demography.add_mass_migration(
        time=admix_gen,
        source="P3",
        dest="P5",
        proportion=alpha
    )

    # Split P3 & P1 from ancestral P13 at time=XA_split_gen
    demography.add_population_split(
        time=XA_split_gen,
        derived=["P3", "P1"],
        ancestral="P13"
    )

    # Split P13 & P5 from ancestral P513 at time=split_gen
    demography.add_population_split(
        time=split_gen,
        derived=["P13", "P5"],
        ancestral="P513"
    )

    # 6.2 Run ancestry simulation
    replicates = msprime.sim_ancestry(
        samples={0: samples, 1: samples, 2: samples, 3: samples, 4: samples},
        sequence_length=L,
        demography=demography,
        recombination_rate=ro,
        num_replicates=num_reps,
        random_seed=seed,
        model=[
            msprime.DiscreteTimeWrightFisher(duration=20),
            msprime.StandardCoalescent()
        ]
    )

    # 6.3 Collect f2 & f3 across replicates
    f2_p3_p1_vals, f2_p3_p5_vals, f2_p1_p5_vals, f3_vals = [], [], [], []
    for ts in replicates:
        # branch‐mode f2 statistics
        a = ts.f2(sample_sets=[ts.samples(3), ts.samples(0)], mode="branch")  # P3–P1
        b = ts.f2(sample_sets=[ts.samples(3), ts.samples(1)], mode="branch")  # P3–P5
        c = ts.f2(sample_sets=[ts.samples(0), ts.samples(1)], mode="branch")  # P1–P5

        # branch‐mode f3: f3(P3; P1, P5)
        #   = ½ [f2(P3,P1) + f2(P3,P5) − f2(P1,P5)]
        f3_val = 0.5 * (a + b - c)

        # accumulate
        f2_p3_p1_vals.append(a)
        f2_p3_p5_vals.append(b)
        f2_p1_p5_vals.append(c)
        f3_vals.append(f3_val)

    # 6.4 Compute means & SDs, store in results
    results["alpha"].append(alpha)
    results["f2_p3_p1"].append(np.mean(f2_p3_p1_vals))
    results["sd_f2_p3_p1"].append(np.std(f2_p3_p1_vals))
    results["f2_p3_p5"].append(np.mean(f2_p3_p5_vals))
    results["sd_f2_p3_p5"].append(np.std(f2_p3_p5_vals))
    results["f2_p1_p5"].append(np.mean(f2_p1_p5_vals))
    results["sd_f2_p1_p5"].append(np.std(f2_p1_p5_vals))
    results["f3_branch"].append(np.mean(f3_vals))
    results["sd_f3_branch"].append(np.std(f3_vals))

# -----------------------------------------------------------------------------
# 7. Convert to pandas DataFrame & display
# -----------------------------------------------------------------------------
df = pd.DataFrame(results)
print(df.head(10))

# -----------------------------------------------------------------------------
# 8. Save results to disk
# -----------------------------------------------------------------------------
with open("msp_f3_variable_alpha_dict.pkl", "wb") as f:
    pickle.dump(results, f)

df.to_pickle("msp_f3_variable_alpha_df.pkl")
print("Data and pickles saved successfully.")

# -----------------------------------------------------------------------------
# 9. Plot f3 & f2 statistics vs. α
# -----------------------------------------------------------------------------
sns.set(style="whitegrid", font_scale=1.2)
palette = sns.color_palette("colorblind", 4)

fig, ax = plt.subplots(figsize=(10, 6))
# Bar: f3 with error bars
ax.bar(
    df["alpha"], df["f3_branch"],
    width=0.02, yerr=df["sd_f3_branch"],
    color=palette[0], capsize=4,
    label="f3(P3; P1, P5)"
)
# Lines: f2s with error bars
ax.errorbar(df["alpha"], df["f2_p3_p1"], yerr=df["sd_f2_p3_p1"],
            fmt="-o", color=palette[1], capsize=3, label="f2(P3,P1)")
ax.errorbar(df["alpha"], df["f2_p3_p5"], yerr=df["sd_f2_p3_p5"],
            fmt="-s", color=palette[2], capsize=3, label="f2(P3,P5)")
ax.errorbar(df["alpha"], df["f2_p1_p5"], yerr=df["sd_f2_p1_p5"],
            fmt="-d", color=palette[3], capsize=3, label="f2(P1,P5)")

ax.set_xlabel("Admixture proportion α (P5 → P3)")
ax.set_ylabel("Branch‐based statistic")
ax.set_title("f3 and f2 statistics vs. admixture proportion")
ax.invert_xaxis()                # show 0.5→0.0 left→right
ax.legend(loc="upper right")
plt.tight_layout()

# Save & show figure
fig.savefig(os.path.join(outdir, "f3_f2_vs_alpha.png"), dpi=300)
plt.show()
