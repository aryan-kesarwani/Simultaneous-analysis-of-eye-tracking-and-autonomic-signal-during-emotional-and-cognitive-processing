
# ===================================================================
# Eye Tracking Analysis Pipeline - Complete Implementation
# Python equivalent of MATLAB getValues_sj.m function
# For Late-Life Depression and Emotional Face-Word Stroop Research
# ===================================================================

import numpy as np
import pandas as pd
import scipy.io
from scipy import interpolate
from scipy.ndimage import uniform_filter1d
import os
import glob
import warnings
import matplotlib.pyplot as plt
import seaborn as sns

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# STEP 1: Configuration and Setup
# ==============================


def struct_to_dataframe(struct_array):
    """
    Convert numpy structured array (MATLAB struct) to pandas DataFrame.
    """
    if struct_array is None:
        return None
    if not hasattr(struct_array, 'dtype') or struct_array.dtype.names is None:
        return None
    data_dict = {}
    for name in struct_array.dtype.names:
        # Extract field and flatten
        col_data = struct_array[name]
        # MATLAB may store data as arrays inside fields, convert to flat lists if possible
        try:
            # If array of size 1, extract scalars
            if isinstance(col_data, np.ndarray):
                if col_data.dtype == 'O':
                    # If cell arrays, convert elements
                    col_data = [elem.item() if hasattr(elem, 'item') else elem for elem in col_data]
                else:
                    col_data = col_data.flatten()
        except Exception:
            pass
        data_dict[name] = col_data
    return pd.DataFrame(data_dict)



def get_dv(subj_dir, dv, tw=None):
    """
    Extract dependent variable data from subject folder.
    Supports beh.mat, info.mat, pupil data, saccade matrices.
    Implements time window averaging similarly to MATLAB getDV.m.

    Parameters
    ----------
    subj_dir : str
        Path to subject folder
    dv : str
        Dependent variable key (e.g., 'RT', 'Pupil_l', 'Rate', etc.)
    tw : tuple or None
        Time window as (start, end); if None, return full data

    Returns
    -------
    dVar : np.ndarray
        Data matrix for the dependent variable
    sub_info : dict or None
        Subject info loaded from info.mat (may be needed for epoch handling)
    """

    print(subj_dir)
    beh_files = glob.glob(os.path.join(subj_dir, '*beh.mat'))
    if not beh_files:
        raise FileNotFoundError(f"No beh.mat file found in {subj_dir}")
    beh_data = scipy.io.loadmat(beh_files[0], struct_as_record=False, squeeze_me=True)
    if 'behData' in beh_data:
        beh_struct = beh_data['behData']
        beh_df = struct_to_dataframe(beh_struct)
        if beh_df is None or beh_df.empty:
            beh_df = None
    else:
        beh_df = None
        print("No behData found in beh.mat")


    info_files = glob.glob(os.path.join(subj_dir, '*info.mat'))
    if not info_files:
        raise FileNotFoundError(f"No info.mat file found in {subj_dir}")
    info_data = scipy.io.loadmat(info_files[0], struct_as_record=False, squeeze_me=True)
    sub_info = info_data.get('subInfo', None)

    if beh_df is not None and dv in beh_df.columns:
        dVar = beh_df[dv].values
    elif sub_info is not None and hasattr(sub_info, dv):
        dVar = getattr(sub_info, dv)
    elif dv == 'Pupil_l':
        files = glob.glob(os.path.join(subj_dir, '*l_epbcData.mat'))
        if not files:
            raise FileNotFoundError(f"No left eye baseline corrected data found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['l_epbcData']
    elif dv == 'Pupil_r':
        files = glob.glob(os.path.join(subj_dir, '*r_epbcData.mat'))
        if not files:
            raise FileNotFoundError(f"No right eye baseline corrected data found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['r_epbcData']
    elif dv == 'tonicPupil_l':
        files = glob.glob(os.path.join(subj_dir, '*l_epData.mat'))
        if not files:
            raise FileNotFoundError(f"No left eye raw data found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['l_epData']
    elif dv == 'tonicPupil_r':
        files = glob.glob(os.path.join(subj_dir, '*r_epData.mat'))
        if not files:
            raise FileNotFoundError(f"No right eye raw data found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['r_epData']
    elif dv == 'Rate':
        files = glob.glob(os.path.join(subj_dir, '*RateMatx.mat'))
        if not files:
            raise FileNotFoundError(f"No RateMatx.mat found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['RateMatx']
        if tw is not None:
            dVar = dVar * (tw[1] - tw[0])
    elif dv == 'Amp':
        files = glob.glob(os.path.join(subj_dir, '*AmpMatx.mat'))
        if not files:
            raise FileNotFoundError(f"No AmpMatx.mat found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['AmpMatx']
    elif dv == 'Vel':
        files = glob.glob(os.path.join(subj_dir, '*VelMatx.mat'))
        if not files:
            raise FileNotFoundError(f"No VelMatx.mat found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['VelMatx']
    elif dv in ['Ang', 'Cos', 'Sin']:
        files = glob.glob(os.path.join(subj_dir, '*AngMatx.mat'))
        if not files:
            raise FileNotFoundError(f"No AngMatx.mat found in {subj_dir}")
        data = scipy.io.loadmat(files[0])
        dVar = data['AngMatx']
    else:
        raise ValueError(f'Unsupported DV: {dv}')

    # Apply time window averaging if specified
    if tw is not None and dVar.ndim > 1 and sub_info and hasattr(sub_info, "epoch"):
        epoch_start = sub_info.epoch[0] if hasattr(sub_info.epoch, '__len__') else sub_info.epoch
        tw_adj = (tw[0] - epoch_start, tw[1] - epoch_start)
        tw_start = max(0, int(tw_adj[0]))
        tw_end = min(dVar.shape[0], int(tw_adj[1]))
        if tw_end > tw_start:
            dVar = np.nanmean(dVar[tw_start:tw_end, :], axis=0)
        else:
            warnings.warn("Time window yields empty slice; skipping time-window averaging.")

    return dVar, sub_info


def get_iv(subj_dir, iv=None, iv_val=None):
    """
    Retrieve trial indices and labels based on independent variables and quality flags.

    Parameters
    ----------
    subj_dir : str
        Subject folder path
    iv : list or None
        Independent variable names
    iv_val : list or None
        Manual condition definitions

    Returns
    -------
    labels : pandas.DataFrame
        Condition label table
    rm_idx : list of boolean arrays
        Trial indices arrays after filtering with rmFlag.overall
    """
    # Load rmFlag first so we can proceed even if behData is missing when iv=None
    rm_flag_files = glob.glob(os.path.join(subj_dir, '*rmFlag.mat'))
    if not rm_flag_files:
        raise FileNotFoundError(f"No rmFlag.mat found in {subj_dir}")
    flag_data = scipy.io.loadmat(rm_flag_files[0], struct_as_record=False, squeeze_me=True)
    if 'rmFlag' not in flag_data:
        raise ValueError("rmFlag missing in rmFlag.mat")
    rmFlag_struct = flag_data['rmFlag']
    overall_field = None
    # Handle numpy structured array
    if isinstance(rmFlag_struct, np.ndarray) and getattr(rmFlag_struct, 'dtype', None) is not None and rmFlag_struct.dtype.names is not None:
        if 'overall' in rmFlag_struct.dtype.names:
            overall_field = rmFlag_struct['overall']
    # Handle MATLAB mat_struct object with attribute access
    if overall_field is None and hasattr(rmFlag_struct, 'overall'):
        overall_field = rmFlag_struct.overall
    if overall_field is None:
        warnings.warn("rmFlag.overall missing, including all trials")
        overall_flag = np.zeros(0, dtype=bool)
    else:
        overall_flag = np.array(overall_field).astype(bool).flatten()

    # Load beh.mat and convert if present
    beh_files = glob.glob(os.path.join(subj_dir, '*beh.mat'))
    beh_df = None
    if beh_files:
        beh_data = scipy.io.loadmat(beh_files[0], struct_as_record=False, squeeze_me=True)
        beh_struct = beh_data.get('behData', None)
        beh_df = struct_to_dataframe(beh_struct)

    # Determine number of trials
    num_trials = len(overall_flag) if overall_flag.size > 0 else (0 if beh_df is None else len(beh_df))
    overall_mask = ~overall_flag if overall_flag.size > 0 else np.ones(num_trials, dtype=bool)

    if iv is None or iv == "none" or len(iv) == 0:
        rm_idx = [overall_mask]
        labels = pd.DataFrame()
    else:
        if beh_df is None or beh_df.empty:
            raise ValueError("behData missing or empty in beh.mat; cannot build IV-based conditions. Use iv=None or provide valid behData.")
        if isinstance(iv, str):
            iv = [iv]

        if iv_val is None:
            groupings = beh_df[iv]
            unique_combos = groupings.drop_duplicates().reset_index(drop=True)
            labels = unique_combos
            rm_idx = []
            for _, combo_row in unique_combos.iterrows():
                mask = np.ones(num_trials, dtype=bool)
                for col in iv:
                    mask &= (beh_df[col] == combo_row[col])
                rm_idx.append(overall_mask & mask.values)
        else:
            rm_idx = []
            for val_row in iv_val:
                mask = np.ones(num_trials, dtype=bool)
                for col_i, col_name in enumerate(iv):
                    mask &= (beh_df[col_name] == val_row[col_i])
                rm_idx.append(overall_mask & mask.values)
            labels = pd.DataFrame(iv_val, columns=iv)

    return labels, rm_idx


def get_values_sj(subj_dir, dv, tw=None, iv=None, iv_val=None):
    """
    Main function to extract and aggregate dependent variable data
    by condition, replicating MATLAB getValues_sj.m.

    Parameters
    ----------
    subj_dir : str
        Subject folder path
    dv : str
        Dependent variable name
    tw : tuple or None
        Time window (start, end)
    iv : list or None
        Independent variables
    iv_val : list or None
        Manual condition specification

    Returns
    -------
    pandas.DataFrame
        Results table with trial counts, condition labels, and mean values
    """
    dVar, sub_info = get_dv(subj_dir, dv, tw)
    labels, rm_idx = get_iv(subj_dir, iv, iv_val)

    n_conds = len(rm_idx)
    # Prepare output arrays
    if dVar.ndim == 1:
        mean_mat = np.full((n_conds,), np.nan)
        num_trials = np.zeros(n_conds, dtype=int)
        for i in range(n_conds):
            idx_mask = rm_idx[i]
            if np.any(idx_mask):
                mean_mat[i] = np.nanmean(dVar[idx_mask])
                num_trials[i] = np.sum(idx_mask)
    else:
        mean_mat = np.full((n_conds, dVar.shape[0]), np.nan)
        num_trials = np.zeros(n_conds, dtype=int)
        for i in range(n_conds):
            idx_mask = rm_idx[i]
            if np.any(idx_mask):
                mean_mat[i, :] = np.nanmean(dVar[:, idx_mask], axis=1)
                num_trials[i] = np.sum(idx_mask)

    # Data-specific processing similar to MATLAB smoothing/interpolation:
    if tw is not None:
        mean_val = mean_mat
    else:
        if dv == "Rate":
            # Moving average with window 100 (approximate)
            mean_val = uniform_filter1d(mean_mat, size=100, axis=1) if mean_mat.ndim>1 else uniform_filter1d(mean_mat, size=100)
        elif dv in ["Amp", "Vel", "Ang", "Cos", "Sin"]:
            # Interpolation + smoothing over 500 samples
            mean_val = mean_mat.copy()
            if dv == "Cos":
                mean_mat = np.cos(mean_mat)
            elif dv == "Sin":
                mean_mat = np.sin(mean_mat)
            for j in range(mean_mat.shape[0]):
                d_sub = mean_mat[j, :]
                # Fill start/end NaNs
                if np.isnan(d_sub[0]): d_sub[0] = np.nanmean(d_sub)
                if np.isnan(d_sub[-1]): d_sub[-1] = np.nanmean(d_sub)
                valid = ~np.isnan(d_sub)
                if np.sum(valid) > 1:
                    f = interpolate.interp1d(np.where(valid)[0], d_sub[valid], kind='linear', fill_value='extrapolate')
                    d_sub_interp = f(np.arange(len(d_sub)))
                    mean_val[j, :] = uniform_filter1d(d_sub_interp, size=500)
        else:
            mean_val = mean_mat

    # Build output DataFrame
    out_data = {"n": num_trials}
    if not labels.empty:
        for col in labels.columns:
            out_data[col] = labels[col].values
    if mean_val.ndim == 1 or mean_val.shape[1] == 1:
        out_data[dv] = mean_val.flatten()
    else:
        # Expand columns for time series
        for t in range(mean_val.shape[1]):
            out_data[f"{dv}_t{t}"] = mean_val[:, t]

    return pd.DataFrame(out_data)


def main_analysis():
    """
    Main analysis function - modify paths and parameters as needed
    """

    # MODIFY THESE PATHS TO MATCH YOUR DATA
    base_directory = "/path/to/your/eye_tracking/data"  # Change this!
    output_directory = "/path/to/output/folder"         # Change this!

    # Define subject information
    # MODIFY THESE TO MATCH YOUR SUBJECT IDs
    healthy_subjects = [f"HC_{i:03d}" for i in range(1, 36)]      # HC_001 to HC_035
    depression_subjects = [f"LLD_{i:03d}" for i in range(1, 36)]  # LLD_001 to LLD_035
    all_subject_ids = healthy_subjects + depression_subjects

    # Create group labels
    group_labels = {}
    for subj in healthy_subjects:
        group_labels[subj] = 'healthy'
    for subj in depression_subjects:
        group_labels[subj] = 'depression'

    print(f"Analysis setup: {len(all_subject_ids)} subjects")
    print(f"Groups: {len(healthy_subjects)} healthy, {len(depression_subjects)} depression")

    # Define analysis parameters
    dependent_variables = [
        'RT',           # Response time (behavioral)
        'ACC',          # Accuracy (behavioral)  
        'Pupil_l',      # Left pupil (baseline corrected)
        'Pupil_r',      # Right pupil (baseline corrected)
        'Rate',         # Saccade rate
        'Amp',          # Saccade amplitude
        'Vel',          # Saccade velocity
    ]

    # Define experimental conditions
    # MODIFY THESE TO MATCH YOUR EXPERIMENTAL DESIGN
    independent_variables = ['face_emotion', 'word_emotion']

    # Define time windows for analysis (in samples or ms - adjust based on your epoch timing)
    time_windows = {
        'Pupil_l': (0, 2000),   # 0-2000ms post-stimulus
        'Pupil_r': (0, 2000),   # 0-2000ms post-stimulus
        'RT': None,             # No time window for behavioral measures
        'ACC': None,            # No time window for behavioral measures
        'Rate': None,           # Full epoch for saccade measures
        'Amp': None,
        'Vel': None,
    }

    # STEP 2: Data Processing
    # ======================

    print("\nStarting data processing...")

    # Process all subjects and dependent variables
    try:
        all_data = process_all_subjects(
            base_dir=base_directory,
            subject_ids=all_subject_ids,
            dependent_vars=dependent_variables,
            independent_vars=independent_variables,
            iv_values=None,  # Auto-extract all combinations
            time_windows=time_windows,
            group_labels=group_labels
        )

        print(f"\n✓ Data processing complete!")
        print(f"Successfully processed {len(all_data)} dependent variables")

    except Exception as e:
        print(f"❌ Error during data processing: {str(e)}")
        return None

    # STEP 3: Generate Summary Statistics
    # ==================================

    print("\nGenerating summary statistics...")
    summary_stats = create_analysis_summary(all_data, group_var='group')

    print("\nSummary Statistics:")
    print("=" * 50)
    print(summary_stats.to_string(index=False))

    # Save summary to CSV
    summary_path = os.path.join(output_directory, "analysis_summary.csv")
    summary_stats.to_csv(summary_path, index=False)
    print(f"\n✓ Summary saved to: {summary_path}")

    # STEP 4: Create Visualizations
    # ============================

    print("\nCreating visualizations...")

    # Create output directory for plots
    plots_dir = os.path.join(output_directory, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    for dv_name in dependent_variables:
        if dv_name in all_data:
            try:
                print(f"  Creating plots for {dv_name}...")

                # Prepare data for visualization
                viz_data = prepare_visualization_data(all_data, dv_name, 'group')

                # 1. Group comparison plots
                fig1 = plot_group_comparison(viz_data, dv_name, 'group', 'condition')
                fig1.savefig(os.path.join(plots_dir, f"{dv_name}_group_comparison.png"), 
                           dpi=300, bbox_inches='tight')
                plt.close(fig1)

                # 2. Congruency effects (for appropriate measures)
                if dv_name in ['RT', 'ACC']:  # Behavioral measures
                    fig2 = plot_congruency_effects(viz_data, dv_name, 'group')
                    fig2.savefig(os.path.join(plots_dir, f"{dv_name}_congruency_effects.png"), 
                               dpi=300, bbox_inches='tight')
                    plt.close(fig2)

                # 3. Time series plots (for pupil and saccade measures)
                if dv_name in ['Pupil_l', 'Pupil_r', 'Rate', 'Amp', 'Vel']:
                    fig3 = plot_time_series(all_data, dv_name, 'group', 'condition')
                    if fig3 is not None:
                        fig3.savefig(os.path.join(plots_dir, f"{dv_name}_time_series.png"), 
                                   dpi=300, bbox_inches='tight')
                        plt.close(fig3)

            except Exception as e:
                print(f"    ❌ Error creating plots for {dv_name}: {str(e)}")

    print(f"\n✓ Visualizations saved to: {plots_dir}")

    # STEP 5: Statistical Analysis Preparation
    # =======================================

    print("\nPreparing data for statistical analysis...")

    # Save individual DV datasets for further analysis
    data_dir = os.path.join(output_directory, "processed_data")
    os.makedirs(data_dir, exist_ok=True)

    for dv_name, data in all_data.items():
        data_path = os.path.join(data_dir, f"{dv_name}_data.csv")
        data.to_csv(data_path, index=False)
        print(f"  Saved {dv_name} data: {data_path}")

    # Create combined dataset for mixed-effects modeling
    combined_data_list = []
    for dv_name, data in all_data.items():
        data_long = data.copy()
        data_long['dependent_variable'] = dv_name

        # Handle array data by taking mean or specific measures
        if hasattr(data_long[dv_name].iloc[0], '__len__') and not isinstance(data_long[dv_name].iloc[0], str):
            # For time series data, compute summary statistics
            data_long[f'{dv_name}_mean'] = data_long[dv_name].apply(lambda x: np.nanmean(x) if hasattr(x, '__len__') else x)
            data_long[f'{dv_name}_peak'] = data_long[dv_name].apply(lambda x: np.nanmax(x) if hasattr(x, '__len__') else x)
            data_long[f'{dv_name}_std'] = data_long[dv_name].apply(lambda x: np.nanstd(x) if hasattr(x, '__len__') else x)

        combined_data_list.append(data_long)

    # Save combined dataset
    combined_path = os.path.join(data_dir, "combined_data.csv")
    if combined_data_list:
        combined_df = pd.concat(combined_data_list, ignore_index=True)
        combined_df.to_csv(combined_path, index=False)
        print(f"  Saved combined data: {combined_path}")

    print("\n" + "="*60)
    print("✓ ANALYSIS COMPLETE!")
    print("="*60)
    print(f"Results saved to: {output_directory}")
    print("\nNext steps:")
    print("1. Review summary statistics and plots")
    print("2. Run statistical tests (t-tests, ANOVA, mixed-effects models)")
    print("3. Examine effect sizes and clinical significance")
    print("4. Consider machine learning classification approaches")

    return all_data, summary_stats

# STEP 6: Usage Instructions
# =========================

if __name__ == "__main__":
    # BEFORE RUNNING:
    # 1. Install required packages: pip install numpy pandas scipy matplotlib seaborn
    # 2. Set correct paths in main_analysis() function
    # 3. Verify your .mat file structure matches the expected format
    # 4. Adjust subject IDs and experimental variables to match your data

    print("Eye Tracking Analysis Pipeline")
    print("==============================")
    print("Make sure to modify paths and parameters before running!")
    print("\nTo run analysis:")
    print("1. Update base_directory and output_directory")
    print("2. Adjust subject_ids to match your data")
    print("3. Verify independent_variables match your experiment")
    print("4. Run: python this_script.py")

    # Uncomment the next line to run the analysis
    # results = main_analysis()
