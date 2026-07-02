# Stage 5: GROMACS MD Setup
# Paste into Python Script node 5, connected after ADMET Ranking
# Requires GROMACS and ACPYPE installed
# Update paths before running

import knime.scripting.io as knio
import pandas as pd
import subprocess
import os
from rdkit import Chem

WORK_DIR = "/Users/sephy/EGFR/md_runs"
ACPYPE   = "/opt/anaconda3/envs/egfr_knime/bin/acpype"
os.makedirs(WORK_DIR, exist_ok=True)

MDP_EM = """
integrator    = steep
emtol         = 1000.0
emstep        = 0.01
nsteps        = 50000
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
pbc           = xyz
"""

input_df = knio.input_tables[0].to_pandas()
top_hits = input_df.head(3).reset_index(drop=True)
results  = []

for idx, row in top_hits.iterrows():
    cdir = os.path.join(WORK_DIR, f"compound_{idx}")
    os.makedirs(cdir, exist_ok=True)

    with open(os.path.join(cdir, "em.mdp"), "w") as f:
        f.write(MDP_EM)

    docked_sdf = row.get('docked_sdf')
    if not docked_sdf or not os.path.exists(str(docked_sdf)):
        results.append({**row.to_dict(), 'md_dir': cdir, 'md_status': 'no_sdf'})
        continue

    mol = Chem.SDMolSupplier(docked_sdf, removeHs=False)[0]
    if mol is None:
        results.append({**row.to_dict(), 'md_dir': cdir, 'md_status': 'mol_read_failed'})
        continue

    lig_pdb = os.path.join(cdir, "ligand.pdb")
    Chem.MolToPDBFile(mol, lig_pdb)

    acpype_cmd = [
        ACPYPE,
        "-i", lig_pdb,
        "-c", "gas",
        "-a", "gaff2",
        "-o", "gmx",
        "-n", "0"
    ]

    r = subprocess.run(acpype_cmd, cwd=cdir, capture_output=True, text=True, timeout=120)

    if r.returncode not in [0, 19]:
        print(f"Compound {idx} ACPYPE failed: {r.stderr[:200]}")
        results.append({**row.to_dict(), 'md_dir': cdir, 'md_status': 'acpype_failed'})
        continue

    results.append({**row.to_dict(), 'md_dir': cdir, 'md_status': 'topology_ready'})
    print(f"Compound {idx}: topology ready at {cdir}")

output_df = pd.DataFrame(results)
knio.output_tables[0] = knio.Table.from_pandas(output_df)
print(f"GROMACS setup complete for {len(output_df)} compounds")
