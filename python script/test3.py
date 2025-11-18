from eye_analysis import get_values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

project_dir = "/Users/aryansmac/Documents/PROJECT/TEEP/TEEP'25/LLD_CTRL_MAIN/Pupil-size-and-microsaccade-analysis-main/new_data/mix_data"

df = get_values(project_dir, dv="Pupil_l", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df2 = get_values(project_dir, dv="Pupil_r", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df3 = get_values(project_dir, dv="Rate", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df4 = get_values(project_dir, dv="Amp", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df5 = get_values(project_dir, dv="Vel", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df6 = get_values(project_dir, dv="Ang", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df7 = get_values(project_dir, dv="Cos", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df8 = get_values(project_dir, dv="Sin", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)

# Use LONG (tidy) format: one 'Value' column and a 'DV' label
def to_long(dfin, dv_name):
    keep = [c for c in dfin.columns if c not in {"Pupil_l","Pupil_r","Rate","Amp","Vel","Ang","Cos","Sin"}]
    return dfin[keep + [dv_name]].rename(columns={dv_name: "Value"}).assign(DV=dv_name)

long_df = pd.concat([
    to_long(df,  "Pupil_l"),
    to_long(df2, "Pupil_r"),
    to_long(df3, "Rate"),
    to_long(df4, "Amp"),
    to_long(df5, "Vel"),
    to_long(df6, "Ang"),
    to_long(df7, "Cos"),
    to_long(df8, "Sin"),
], ignore_index=True)

# explode each timecourse into one row per time index
long_df = long_df.assign(ValueSeq=long_df["Value"].apply(lambda x: list(np.asarray(x).ravel())))
long_time = long_df.explode("ValueSeq", ignore_index=True).rename(columns={"ValueSeq": "Value"})
long_time["t_idx"] = long_time.groupby(["Subj", "DV"]).cumcount()

print(long_time.head(20))