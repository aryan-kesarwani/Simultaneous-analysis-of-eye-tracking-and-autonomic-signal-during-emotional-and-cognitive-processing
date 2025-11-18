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

# Global plotting style and quality
sns.set_theme(style="whitegrid")
sns.set_context("talk", font_scale=1.2)
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.dpi"] = 150
plt.rcParams["axes.linewidth"] = 1.0

# ----- Helpers -----
def to_array(x):
    if isinstance(x, np.ndarray):
        return x.astype(float)
    if isinstance(x, (list, tuple)):
        return np.asarray(x, dtype=float)
    return np.asarray([x], dtype=float)

def analyze_and_plot(df_in: pd.DataFrame, dv_name: str):
    if df_in is None or df_in.empty or dv_name not in df_in.columns:
        return
    df = df_in.copy()
    # Ensure arrays
    df["series"] = df[dv_name].apply(to_array)
    # Trim to common min length for consistent aggregation
    min_len = int(min(s.size for s in df["series"])) if len(df) else 0
    if min_len == 0:
        return
    df["series"] = df["series"].apply(lambda s: s[:min_len])
    df["group"] = df["group"].astype(int)

    # Aggregate by group over time (mean ± SEM)
    group_means = {}
    group_sems = {}
    for g, gdf in df.groupby("group"):
        mat = np.vstack(gdf["series"].values)
        group_means[g] = np.nanmean(mat, axis=0)
        denom = np.sqrt(np.maximum(np.sum(~np.isnan(mat), axis=0), 1))
        group_sems[g] = np.nanstd(mat, axis=0, ddof=1) / denom

    t = np.arange(min_len)
    # Timecourse: mean ± SEM per group
    plt.figure(figsize=(12, 8))
    colors = {1: (0.0, 0.5, 0.8), 2: (0.8, 0.1, 0.0)}
    for g in sorted(group_means.keys()):
        m = group_means[g]
        s = group_sems[g]
        c = colors.get(g, (0.3, 0.3, 0.3))
        plt.plot(t, m, label=f"Group {g}", color=c, linewidth=2.5)
        plt.fill_between(t, m - s, m + s, color=c, alpha=0.2, linewidth=0)
    plt.xlabel("Time (samples)")
    plt.ylabel(dv_name)
    plt.title(f"{dv_name} dynamics by group")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"analysis_{dv_name}_timecourse.png", dpi=400, bbox_inches="tight")
    plt.close()

    # Per-subject mean over time: boxplot by group
    per_subject = (
        pd.DataFrame({
            "Subj": df["Subj"].values,
            "group": df["group"].values,
            dv_name: [np.nanmean(s) for s in df["series"].values],
        })
    )
    plt.figure(figsize=(7, 6))
    sns.boxplot(data=per_subject, x="group", y=dv_name, linewidth=1.2)
    sns.stripplot(data=per_subject, x="group", y=dv_name, color="k", size=4, alpha=0.35)
    plt.title(f"{dv_name}: per-subject mean by group")
    plt.tight_layout()
    plt.savefig(f"analysis_{dv_name}_boxplot.png", dpi=400, bbox_inches="tight")
    plt.close()

# ----- Run per-DV analysis (no concatenation) -----
name_to_df = [
    ("Pupil_l", df),
    ("Pupil_r", df2),
    ("Rate",    df3),
    ("Amp",     df4),
    ("Vel",     df5),
    ("Ang",     df6),
    ("Cos",     df7),
    ("Sin",     df8),
]

for dv_name, dfin in name_to_df:
    analyze_and_plot(dfin, dv_name)

print("Saved per-DV plots: analysis_<DV>_timecourse.png and analysis_<DV>_boxplot.png")