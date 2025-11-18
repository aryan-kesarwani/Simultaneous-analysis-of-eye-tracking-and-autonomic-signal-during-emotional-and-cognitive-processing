from eye_analysis import get_values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

project_dir = "/Users/aryansmac/Documents/PROJECT/TEEP/TEEP'25/LLD_CTRL_MAIN/Pupil-size-and-microsaccade-analysis-main/new_data/mix_data"

# Get per-subject timecourse (no time window); group labels from info
df = get_values(project_dir, dv="Pupil_l", iv=["group"], iv_val=None)

# Convert Pupil_l column (list/array) to scalar float
def to_scalar(x):
    if x is None:
        return np.nan
    if isinstance(x, (list, tuple)):
        return float(x[0]) if len(x) else np.nan
    if isinstance(x, np.ndarray):
        return float(x.ravel()[0]) if x.size else np.nan
    return float(x)

# Prepare per-group timecourses
def to_array(x):
    if isinstance(x, np.ndarray):
        return x.astype(float)
    if isinstance(x, (list, tuple)):
        return np.asarray(x, dtype=float)
    # if scalar due to TW, expand to 1-length
    return np.asarray([x], dtype=float)

df["series"] = df["Pupil_l"].apply(to_array)
df["group"] = df["group"].astype(int)

# Ensure all series have the same length by trimming to min length
min_len = int(min(s.size for s in df["series"]))
df["series"] = df["series"].apply(lambda s: s[:min_len])

# Aggregate: mean and SEM across subjects within each group at each timepoint
group_means = {}
group_sems = {}
for g, gdf in df.groupby("group"):
    mat = np.vstack(gdf["series"].values)  # shape: (n_subjects, T)
    group_means[g] = np.nanmean(mat, axis=0)
    denom = np.sqrt(np.maximum(np.sum(~np.isnan(mat), axis=0), 1))
    group_sems[g] = np.nanstd(mat, axis=0, ddof=1) / denom

# Time axis (sample index). If you need ms, scale using sampleRate
t = np.arange(min_len)

# Colors as specified: [[0 .5 .8]; [.8 .1 0]]
col = {
    1: (0.0, 0.5, 0.8),
    2: (0.8, 0.1, 0.0),
}

plt.figure(figsize=(9, 4))
for g in sorted(group_means.keys()):
    m = group_means[g]
    s = group_sems[g]
    c = col.get(g, (0.3, 0.3, 0.3))
    plt.plot(t, m, label=f"Group {g}", color=c, linewidth=2)
    plt.fill_between(t, m - s, m + s, color=c, alpha=0.2, linewidth=0)

plt.xlabel("Time from target onset (ms)")
plt.ylabel("Pupil size (mm)")
plt.title("Pupil_l dynamics by group")
plt.legend()
plt.tight_layout()
plt.savefig("pupil_l_group_timecourse.png", dpi=150)
plt.show()
