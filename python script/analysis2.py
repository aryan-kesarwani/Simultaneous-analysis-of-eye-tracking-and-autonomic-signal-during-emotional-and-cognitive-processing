from eye_analysis import get_values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

project_dir = "/Users/aryansmac/Documents/PROJECT/TEEP/TEEP'25/LLD_CTRL_MAIN/Pupil-size-and-microsaccade-analysis-main/new_data/mix_data"

df = get_values(project_dir, dv="Pupil_l", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df2 = get_values(project_dir, dv="Pupil_r", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df3 = get_values(project_dir, dv="Rate", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df4 = get_values(project_dir, dv="Amp", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df5 = get_values(project_dir, dv="Vel", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df6 = get_values(project_dir, dv="Ang", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df7 = get_values(project_dir, dv="Cos", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)
df8 = get_values(project_dir, dv="Sin", iv=["group", "age", "gender", "GDS_15_total_score"], iv_val=None)

df3["Rate_mean"] = df3["Rate"].apply(
    lambda x: np.nanmean(x) if isinstance(x, (list, np.ndarray)) else x
)

plt.figure(figsize=(8, 6))
sns.scatterplot(x="GDS_15_total_score", y="Rate_mean", hue="group", data=df3)
sns.regplot(x="GDS_15_total_score", y="Rate_mean", data=df3, scatter=False, color="black")
plt.title("Microsaccade Rate vs Depression Score")
plt.savefig(f"analysis_Correlation_with_Depression_Score_(GDS_15).png", dpi=400, bbox_inches="tight")
plt.show()
