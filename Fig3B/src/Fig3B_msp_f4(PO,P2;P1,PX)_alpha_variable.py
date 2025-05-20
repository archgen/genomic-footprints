# F4 Statistic Power Simulations
# Fig 3B
# This script simulates ancestry under a predefined demographic model using msprime and stdpopsim to 
# compute F2 and F4 statistics across a range of admixture proportions, then visualizes how those statistics vary with admixture.

# Williams & Huber 2025 – Genomic Footprints 

import os
import msprime           # for coalescent simulations
import tskit             # for tree sequence statistics
import stdpopsim         # to load realistic demographic models
import numpy as np       # numerical computing
import pandas as pd      # tabular data handling
import matplotlib.pyplot as plt
import seaborn as sns    # enhanced plotting
import demesdraw         # (optional) drawing demographies

# 1. Prepare output directory
outdir = "/Users/mkw5910/Documents/admixture_aDNA_review/figures/f3_X_A_B_model/"
os.makedirs(outdir, exist_ok=True)    # create if not already present

# 2. Load genome and demographic model via stdpopsim
engine  = stdpopsim.get_engine("msprime")            # simulation backend
species = stdpopsim.get_species("HomSap")            # Homo sapiens species
demog   = species.get_demographic_model("AncientEurasia_9K19")
contig  = species.get_contig("chr21", mutation_rate=demog.mutation_rate)
L       = contig.recombination_map.position[-1]      # sequence length
ro      = contig.recombination_map.rate[0]           # recombination rate

# 3. Utility function: coalescence probability up to time t
def probability_of_coalescence(Ne, t):
    """
    Returns probability that two lineages coalesce before generation t,
    given effective population size Ne.
    """
    return 1 - np.exp(-t / (2 * Ne))

# 4. Simulation parameters
split_gen      = 100     # generation at final population split
ne             = 1000    # effective population size
time           = split_gen - 40  # intermediate split time
samples        = 20      # number of samples per population
seed           = 7       # RNG seed for reproducibility
num_replicates = 1       # replicates per admixture proportion
admix_gen      = 5       # generation of admixture event

# 5. Storage for summary statistics
admix_values = np.arange(0, 1.05, 0.05)  # α from 0.0 to 1.0
results = {
    "alpha":             [],
    "f2_PO_P1":          [], "sd_PO_P1": [],
    "f2_P2_PX":          [], "sd_P2_PX": [],
    "f2_PO_PX":          [], "sd_PO_PX": [],
    "f2_P1_P2":          [], "sd_P1_P2": [],
    "f4_branch":         [], "sd_f4":   [],
    "f4_from_f2_branch": []
}

# 6. Loop over admixture proportions
for alpha in admix_values:
    # 6.1 Build demographic model in msprime
    demography = msprime.Demography()
    demography.add_population(name="P1",  initial_size=ne)   # pop0
    demography.add_population(name="P2",  initial_size=ne)   # pop1
    demography.add_population(name="P1P2",initial_size=ne)   # pop2
    demography.add_population(name="PX",  initial_size=ne)   # pop3 (admixed)
    demography.add_population(name="PO",  initial_size=ne)   # pop4 (outgroup)
    demography.add_population(name="PANC",initial_size=ne)   # pop5 (ancestral)
    # Admixture event at admix_gen: PX ← α·P1 + (1−α)·P2
    demography.add_admixture(
        time=admix_gen,
        derived="PX",
        ancestral=["P1", "P2"],
        proportions=[alpha, 1 - alpha]
    )
    # Split P1/P2 from P1P2 at time 'time'
    demography.add_population_split(time=time, derived=["P1","P2"], ancestral="P1P2")
    # Split P1P2/PO from PANC at split_gen
    demography.add_population_split(time=split_gen, derived=["P1P2","PO"], ancestral="PANC")

    # 6.2 Run ancestry simulations
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

    # 6.3 Compute F2/F4 per replicate
    f2_PO_P1, f2_P2_PX, f2_PO_PX, f2_P1_P2 = [], [], [], []
    f4_direct, f4_from_f2 = [], []

    for ts in replicates:
        # F4 statistic via branch mode: F4(PO, P2; P1, PX)
        f4_val = ts.f4(
            sample_sets=[ts.samples(4), ts.samples(1), ts.samples(0), ts.samples(3)],
            mode="branch"
        )
        # Four F2 pairs
        a = ts.f2(sample_sets=[ts.samples(4), ts.samples(0)], mode="branch")  # PO-P1
        b = ts.f2(sample_sets=[ts.samples(1), ts.samples(3)], mode="branch")  # P2-PX
        c = ts.f2(sample_sets=[ts.samples(4), ts.samples(3)], mode="branch")  # PO-PX
        d = ts.f2(sample_sets=[ts.samples(0), ts.samples(1)], mode="branch")  # P1-P2

        # F4 from F2 relationships: ½(c + d − a − b)
        f4_via_f2 = 0.5 * (c + d - a - b)

        # Accumulate per‐replicate values
        f4_direct.append(f4_val)
        f4_from_f2.append(f4_via_f2)
        f2_PO_P1.append(a)
        f2_P2_PX.append(b)
        f2_PO_PX.append(c)
        f2_P1_P2.append(d)

    # 6.4 Aggregate means and standard deviations
    results["alpha"].append(alpha)
    results["f4_branch"].append(np.mean(f4_direct))
    results["sd_f4"].append(np.std(f4_direct))
    results["f4_from_f2_branch"].append(np.mean(f4_from_f2))

    results["f2_PO_P1"].append(np.mean(f2_PO_P1))
    results["sd_PO_P1"].append(np.std(f2_PO_P1))
    results["f2_P2_PX"].append(np.mean(f2_P2_PX))
    results["sd_P2_PX"].append(np.std(f2_P2_PX))
    results["f2_PO_PX"].append(np.mean(f2_PO_PX))
    results["sd_PO_PX"].append(np.std(f2_PO_PX))
    results["f2_P1_P2"].append(np.mean(f2_P1_P2))
    results["sd_P1_P2"].append(np.std(f2_P1_P2))

# 7. Create a DataFrame for easy viewing
df = pd.DataFrame(results)
print(df)

# 8. Plot F2 and F4 as functions of admixture proportion α
sns.set(style="whitegrid", font_scale=1.2)
palette = sns.color_palette("colorblind", 5)

fig, ax = plt.subplots(figsize=(10, 6))
# Bar plot for F4 (direct)
ax.bar(df["alpha"], df["f4_branch"], width=0.02, yerr=df["sd_f4"],
       color=palette[0], capsize=4, label="F4 (branch)")
# Line plots for each F2 pair
ax.errorbar(df["alpha"], df["f2_PO_P1"], yerr=df["sd_PO_P1"],
            fmt="-o", color=palette[1], label="F2(PO,P1)")
ax.errorbar(df["alpha"], df["f2_P2_PX"], yerr=df["sd_P2_PX"],
            fmt="-o", color=palette[2], label="F2(P2,PX)")
ax.errorbar(df["alpha"], df["f2_PO_PX"], yerr=df["sd_PO_PX"],
            fmt="-o", color=palette[3], label="F2(PO,PX)")
ax.errorbar(df["alpha"], df["f2_P1_P2"], yerr=df["sd_P1_P2"],
            fmt="-o", color=palette[4], label="F2(P1,P2)")

ax.set_xlabel("Admixture proportion α (P1 → PX)")
ax.set_ylabel("Branch-based statistic")
ax.set_title("F2 and F4 statistics vs. admixture proportion")
ax.legend(loc="upper right")
ax.invert_xaxis()       # if desired: show 1→0 on the x-axis
plt.tight_layout()
plt.show()
