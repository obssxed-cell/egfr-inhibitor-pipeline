# Stage 2: DeepChem Activity Prediction
# Paste into Python Script node 2, connected after RDKit Prep
# Trains on EGFR ChEMBL IC50 data and predicts activity for all compounds
# Update EGFR_CSV path before running

import knime.scripting.io as knio
import pandas as pd
import deepchem as dc
import numpy as np
import os
os.environ['DGLBACKEND'] = 'pytorch'

EGFR_CSV = "/Users/sephy/EGFR/data/DOWNLOAD-gU8RPQ5Wut7KaKJdHzr2fUYYJcpIjb0ClUND2cUakNk_eq_.csv"

input_df = knio.input_tables[0].to_pandas()
smiles   = input_df['smiles'].tolist()

egfr_df = pd.read_csv(EGFR_CSV, sep=';')
egfr_df = egfr_df[['Smiles', 'pChEMBL Value']].dropna()
egfr_df.columns = ['smiles', 'pchembl']
egfr_df['pchembl'] = pd.to_numeric(egfr_df['pchembl'], errors='coerce')
egfr_df = egfr_df.dropna()
egfr_df['active'] = (egfr_df['pchembl'] >= 6.0).astype(int)

featurizer = dc.feat.CircularFingerprint(size=1024, radius=2)
X_train = featurizer.featurize(egfr_df['smiles'].tolist())
y_train = egfr_df['active'].values.reshape(-1, 1)
train_ds = dc.data.NumpyDataset(X=X_train, y=y_train)

model = dc.models.MultitaskClassifier(n_tasks=1, n_features=1024, layer_sizes=[512, 256])
model.fit(train_ds, nb_epoch=20)

X = featurizer.featurize(smiles)
pred_ds = dc.data.NumpyDataset(X=X, ids=smiles)
preds = model.predict_on_batch(X)
input_df['predicted_active'] = preds[:, 0, 1] if preds.ndim == 3 else preds.flatten()

output_df = input_df[input_df['predicted_active'] > 0.5].reset_index(drop=True)
knio.output_tables[0] = knio.Table.from_pandas(output_df)
print(f"DeepChem screen: {len(output_df)} predicted actives from {len(input_df)}")
