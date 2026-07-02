# Stage 1: RDKit Ligand Preparation
# Paste into Python Script node 1, connected after CSV Reader
# Input column must be named 'Smiles' (ChEMBL default)

import knime.scripting.io as knio
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

input_df = knio.input_tables[0].to_pandas()
smiles_col = 'Smiles'

records = []
for _, row in input_df.iterrows():
    smi = str(row[smiles_col])
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.FragmentParent(mol)
    mol = rdMolStandardize.ChargeParent(mol)

    mw   = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    if sum([mw > 500, logp > 5, hbd > 5, hba > 10]) > 1:
        continue

    records.append({
        'smiles': Chem.MolToSmiles(mol),
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'HBD': hbd,
        'HBA': hba,
        'TPSA': round(tpsa, 2),
        'RotBonds': rotb,
    })

output_df = pd.DataFrame(records)
knio.output_tables[0] = knio.Table.from_pandas(output_df)
print(f"Lipinski filter: {len(output_df)}/{len(input_df)} passed")
