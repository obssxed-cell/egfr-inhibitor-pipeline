# Stage 3: GNINA Molecular Docking
# Paste into Python Script node 3, connected after DeepChem Screen
# Requires Docker with gnina/gnina image pulled
# Update RECEPTOR and REF_LIG paths before running

import knime.scripting.io as knio
import pandas as pd
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem

RECEPTOR   = "/EGFR/data/4HJO_prepared.pdb"
REF_LIG    = "/EGFR/data/AQ4_ref.pdb"
OUTPUT_DIR = "/Users/sephy/EGFR/gnina_outputs"
DOCKER     = "/usr/local/bin/docker"
os.makedirs(OUTPUT_DIR, exist_ok=True)

input_df = knio.input_tables[0].to_pandas()
input_df = input_df.nlargest(20, 'predicted_active').reset_index(drop=True)

results = []
for idx, row in input_df.iterrows():
    mol = Chem.MolFromSmiles(row['smiles'])
    if mol is None:
        continue

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        continue
    AllChem.MMFFOptimizeMolecule(mol)

    lig_path = os.path.join(OUTPUT_DIR, f"lig_{idx}.sdf")
    out_path = os.path.join(OUTPUT_DIR, f"docked_{idx}.sdf")
    writer = Chem.SDWriter(lig_path)
    writer.write(mol)
    writer.close()

    lig_docker = lig_path.replace("/Users/sephy/EGFR", "/EGFR")
    out_docker = out_path.replace("/Users/sephy/EGFR", "/EGFR")

    cmd = [
        DOCKER, "run", "--rm",
        "-v", "/Users/sephy/EGFR:/EGFR",
        "gnina/gnina", "gnina",
        "-r", RECEPTOR,
        "-l", lig_docker,
        "--autobox_ligand", REF_LIG,
        "--exhaustiveness", "8",
        "--num_modes", "5",
        "--cnn", "fast",
        "-o", out_docker,
        "--cpu", "4",
        "--quiet"
    ]

    cnn_affinity = None
    cnn_score    = None
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if os.path.exists(out_path):
            supp = Chem.SDMolSupplier(out_path, removeHs=False)
            for dm in supp:
                if dm:
                    props = dm.GetPropsAsDict()
                    cnn_affinity = props.get('CNNaffinity')
                    cnn_score    = props.get('CNNscore')
                    break
    except Exception as e:
        print(f"Compound {idx} failed: {e}")

    results.append({
        **row.to_dict(),
        'cnn_affinity': cnn_affinity,
        'cnn_score':    cnn_score,
        'docked_sdf':   out_path,
    })
    print(f"Compound {idx}: CNNaffinity={cnn_affinity}, CNNscore={cnn_score}")

output_df = pd.DataFrame(results)
output_df = output_df.dropna(subset=['cnn_affinity']).reset_index(drop=True)
knio.output_tables[0] = knio.Table.from_pandas(output_df)
print(f"Docking complete: {len(output_df)} compounds docked")
