# EGFR Inhibitor Discovery Pipeline

Automated multi-stage virtual screening pipeline built in KNIME Analytics Platform 
for the computational discovery of potential EGFR (Epidermal Growth Factor Receptor) 
inhibitors.

## Pipeline Overview

**Stage 1 — Ligand Preparation (RDKit)**  
Sanitises and standardises 58,847 ChEMBL EGFR bioactivity records. Applies 
Lipinski's Rule of Five filtering and computes physicochemical descriptors. 
Output: 19,929 drug-like compounds.

**Stage 2 — Activity Prediction (DeepChem)**  
Trains a MultitaskClassifier neural network on the EGFR ChEMBL IC50 data using 
ECFP fingerprints (1024 bits, radius 2). Labels compounds with pChEMBL >= 6 as 
active and predicts activity across the filtered library. Output: 14,846 predicted 
EGFR actives.

**Stage 3 — Molecular Docking (GNINA)**  
Docks top compounds into the EGFR ATP-binding pocket using the 4HJO crystal 
structure (erlotinib-bound, 1.5 Å) as receptor. GNINA's CNN ensemble scoring 
generates CNNaffinity and CNNscore for each pose. Requires Docker.

**Stage 4 — ADMET Flagging**  
Rule-based hERG cardiotoxicity and bioavailability flags. Composite rank score 
combining CNN affinity, predicted activity and ADMET flags.

**Stage 5 — GROMACS MD Setup**  
Generates GAFF2 force field topology files for top hits using ACPYPE. Output 
directories contain ligand.itp, ligand.gro and ligand.top ready for MD simulation.


## Resulta

### Example Output

Top 5 hits from screening 58,847 ChEMBL EGFR compounds:

| Rank | MW | LogP | CNNaffinity | CNNscore | rank_score |
|------|----|------|-------------|----------|------------|
| 1 | 589.29 | 4.40 | 8.665 | 0.297 | 4.633 |
| 2 | 541.27 | 4.39 | 8.376 | 0.455 | 4.588 |
| 3 | 443.18 | 5.44 | 7.797 | 0.835 | 4.298 |
| 4 | 384.20 | 4.95 | 7.381 | 0.601 | 4.190 |
| 5 | 443.18 | 5.44 | 7.464 | 0.704 | 4.132 |

## Requirements

- KNIME Analytics Platform 5.11
- Docker Desktop (for GNINA)
- GROMACS 2026.3
- conda

## Installation

```bash
conda env create -f environment.yml
conda activate egfr_knime
```

Point KNIME at the egfr_knime environment under:  
Preferences → KNIME → Conda → set path to your Anaconda installation  
Preferences → KNIME → Python (legacy) → Python 3 Conda environment → egfr_knime

## Data

Download EGFR (CHEMBL203) bioactivity data from ChEMBL:  
https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL203/

Filter by Standard Type = IC50 and download as CSV with semicolon delimiter.

## Receptor Preparation

Download 4HJO from the RCSB PDB. In ChimeraX:

```
open 4HJO.pdb
delete solvent
delete :AQ4
addh
save 4HJO_prepared.pdb
```

Extract reference ligand:

```
open 4HJO.pdb
select #1/A:1001
save AQ4_ref.pdb selectedOnly true
```

## GNINA

Pull the Docker image:

```bash
docker pull gnina/gnina
```

## Results

Top docking hit: CNNaffinity 8.67 kcal/mol  
Best balanced hit: MW 361, LogP 4.04, CNNaffinity 7.02, CNNscore 0.88

## Next Steps

100ns GROMACS MD simulation on top 3 hits on HPC cluster. MM-PBSA binding 
free energy calculations.

## Target

EGFR (Epidermal Growth Factor Receptor, CHEMBL203) — a receptor tyrosine kinase 
overexpressed in multiple cancers and the target of approved inhibitors including 
erlotinib, gefitinib and osimertinib.

## Author

Instagram: sephyscript
GitHub: obssxed-cell
