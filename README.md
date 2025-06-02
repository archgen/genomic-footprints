# Genomic Footprints

This repository contains simulation scripts used to explore the behavior of f₃ and f₄ statistics under a variety of demographic scenarios. The scripts generate the panels shown in Figures 2 and 3 of our study by running coalescent simulations (via msprime/stdpopsim), computing branch-based f₂/f₃/f₄ statistics, and plotting the results. Figures are from Williams & Huber 2025, The genomic footprints of migration: how ancient DNA is revealing our history of mobility.

Directory structure:

├── Fig2A-B/  
│   └── src/  
│       └── Fig2A-B_msp_f3_(P3;P1,P5)_alpha_variable.py  
├── Fig2C/  
│   └── src/  
│       └── Fig2C_msp_f3_(P3;P1,P5)_alpha_variable.py  
├── Fig2D/  
│   └── src/  
│       └── Fig2D_f3(p3-p1)(p3-p5).py  
├── Fig3B/  
│   └── src/  
│       └── Fig3B_msp_f4(PO,P2;P1,PX)_alpha_variable.py  
└── Fig3C/  
    └── src/  
        └── Fig3C_msp_f4(PO,P2;P1,PX)_splitTime.py  

## Prerequisites

- Python 3.7 or newer  
- The following Python packages (e.g. via pip or conda):

bash pip install msprime tskit stdpopsim demesdraw numpy pandas seaborn matplotlib


## Usage

1. Clone this repository.  
2. Install dependencies as above.  
3. For each figure (2A-B, 2C, 2D, 3B, 3C) cd into its directory and run the script in `src/`.  
4. Output figures and data will be saved to the paths defined in each script (modify the `outdir` variable at the top of a script if needed).

## License
This project is licensed under the [Creative Commons - Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/legalcode "Creative Commons - Attribution 4.0 International License (CC BY 4.0)") - see the [LICENSE](LICENSE) file for details.
     

