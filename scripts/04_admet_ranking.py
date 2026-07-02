# Stage 4: ADMET Flagging and Composite Ranking
# Paste into Python Script node 4, connected after GNINA Docking

import knime.scripting.io as knio
import pandas as pd

input_df = knio.input_tables[0].to_pandas()

input_df['hERG_risk_flag'] = (
    (input_df['LogP'] > 3.5) & (input_df['MW'] > 400)
).astype(int)

input_df['poor_bioavailability_flag'] = (
    (input_df['TPSA'] > 140) | (input_df['RotBonds'] > 10)
).astype(int)

input_df['rank_score'] = (
    input_df['cnn_affinity'] * 0.5 +
    input_df['predicted_active'] * 0.3 +
    (1 - input_df['hERG_risk_flag']) * 0.1 +
    (1 - input_df['poor_bioavailability_flag']) * 0.1
)

output_df = input_df.sort_values('rank_score', ascending=False).reset_index(drop=True)
knio.output_tables[0] = knio.Table.from_pandas(output_df)
print(f"ADMET flagging complete. Top rank score: {output_df['rank_score'].iloc[0]:.3f}")
