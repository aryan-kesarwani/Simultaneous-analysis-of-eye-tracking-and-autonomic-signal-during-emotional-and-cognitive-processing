import os
import numpy as np
import pandas as pd
from scipy.io import loadmat
import h5py
from scipy.interpolate import interp1d


def get_values(project_dir, dv, tw=None, iv=None, iv_val=None, require_beh_info=False, return_meta=False):
    """
    Main function: loops over subjects and collects DV values.
    """
    # Sort by directory name to avoid comparing DirEntry objects directly
    prec_dir = sorted(
        [d for d in os.scandir(os.path.join(project_dir, "derivatives", "prec")) if d.is_dir()],
        key=lambda d: d.name,
    )
    op_sub = []

    for n_subj, subj_entry in enumerate(prec_dir, start=1):
        # Each entry in prec_dir already points to a subject directory like sub-01
        subj_dir = subj_entry.path
        #print("Directory: ", subj_dir)

        # Load info.mat
        info_file = [f for f in os.scandir(subj_dir) if f.name.endswith("info.mat")][0]
        #print("Info file: ", info_file.path)
        sub_info = loadmat(info_file.path, squeeze_me=True, struct_as_record=False)["subInfo"]

        if hasattr(sub_info, "olFlag") and sub_info.olFlag:
            continue

        op = get_values_sj(subj_dir, dv, tw, iv, iv_val, require_beh_info=require_beh_info, return_meta=return_meta)
        if op is None or op.empty:
            continue

        op.insert(0, "Subj", n_subj)  # prepend subject ID
        op_sub.append(op)

    if op_sub:
        return pd.concat(op_sub, ignore_index=True)
    return pd.DataFrame()


def get_values_sj(subj_dir, dv, tw, iv, iv_val, require_beh_info=False, return_meta=False):
    """
    Supporting function: extracts DV, filters IV, computes condition means.
    """
    dvar, beh_meta, info_meta = get_dv(subj_dir, dv, tw, require_beh_info=require_beh_info, return_meta=True)
    labels, rm_idx = get_iv(subj_dir, iv, iv_val)

    mean_mat = []
    counts = []
    for mask in rm_idx:
        valid_idx = np.where(mask)[0]
        if valid_idx.size > 0:
            mean_mat.append(np.nanmean(dvar[:, valid_idx], axis=1))
            counts.append(valid_idx.size)
        else:
            mean_mat.append(np.full(dvar.shape[0], np.nan))
            counts.append(0)

    mean_mat = np.vstack(mean_mat) if mean_mat else np.empty((0, dvar.shape[0]))

    # Post-processing based on DV type
    if tw is not None:
        mean_val = mean_mat
    else:
        if dv == "Rate":
            mean_val = pd.DataFrame([
                pd.Series(row).rolling(100, min_periods=1).mean().values for row in mean_mat
            ]).to_numpy()
        elif dv in {"Amp", "Vel", "Ang", "Cos", "Sin"}:
            mean_val = np.copy(mean_mat)
            if dv in {"Cos", "Sin"}:
                fn = np.cos if dv == "Cos" else np.sin
                mean_mat = fn(mean_mat)

            t = mean_mat.shape[1]
            for i, row in enumerate(mean_mat):
                d_sub = np.copy(row)
                nan_mask = np.isnan(d_sub)
                if nan_mask.all():
                    continue
                d_sub[0] = d_sub[-1] = np.nanmean(d_sub)
                t_sub = np.where(~nan_mask)[0]
                try:
                    f = interp1d(t_sub, d_sub[~nan_mask], kind="linear", fill_value="extrapolate")
                    d_interp = f(np.arange(t))
                    d_smooth = pd.Series(d_interp).rolling(500, min_periods=1).mean().values
                    mean_val[i, :] = d_smooth
                except Exception:
                    pass
        else:
            mean_val = mean_mat

    # Construct output DataFrame
    df = labels.copy()
    df.insert(0, "n", counts)
    df[dv] = list(mean_val)
    if return_meta:
        return df, beh_meta, info_meta
    return df


def get_dv(subj_dir, dv, tw, require_beh_info=False, return_meta=False):
    """
    Extract dependent variable matrix from MAT files.
    """
    def find_by_suffix(suffix):
        files = [f for f in os.scandir(subj_dir) if f.is_file() and f.name.endswith(suffix)]
        if not files:
            raise FileNotFoundError(f"No file ending with {suffix} found in {subj_dir}")
        return files[0].path

    def load_mat_variable(path, var):
        # First try scipy (non-v7.3). If key missing, return first non-meta variable.
        try:
            d = loadmat(path, squeeze_me=True, struct_as_record=False)
            data_keys = [k for k in d.keys() if not k.startswith("__")]
            if var in d:
                return d[var]
            if data_keys:
                return d[data_keys[0]]
        except NotImplementedError:
            # MATLAB v7.3 HDF5
            with h5py.File(path, "r") as f:
                if var in f:
                    obj = f[var]
                    if isinstance(obj, h5py.Dataset):
                        return obj[()]
                    if isinstance(obj, h5py.Group):
                        return {k: (v[()] if isinstance(v, h5py.Dataset) else v) for k, v in obj.items()}
                # best-effort: return first dataset/group
                for k, v in f.items():
                    if isinstance(v, h5py.Dataset):
                        return v[()]
                    if isinstance(v, h5py.Group):
                        return {kk: (vv[()] if isinstance(vv, h5py.Dataset) else vv) for kk, vv in v.items()}
        # If we get here, neither reader produced data
        raise KeyError(f"Variable {var} not found in {path}")

    # Optionally pre-load beh/info if requested
    beh_meta = None
    info_meta = None
    if require_beh_info:
        try:
            beh_meta = load_mat_variable(find_by_suffix("beh.mat"), "behData")
        except Exception:
            pass
        try:
            info_meta = load_mat_variable(find_by_suffix("info.mat"), "subInfo")
        except Exception:
            pass

    # Fast paths that do not require beh/info
    if dv == "Pupil_l":
        dvar = load_mat_variable(find_by_suffix("l_epbcData.mat"), "l_epbcData")
    elif dv == "Pupil_r":
        dvar = load_mat_variable(find_by_suffix("r_epbcData.mat"), "r_epbcData")
    elif dv == "tonicPupil_l":
        dvar = load_mat_variable(find_by_suffix("l_epData.mat"), "l_epData")
    elif dv == "tonicPupil_r":
        dvar = load_mat_variable(find_by_suffix("r_epData.mat"), "r_epData")
    elif dv == "Rate":
        dvar = load_mat_variable(find_by_suffix("RateMatx.mat"), "RateMatx")
        if tw is not None:
            dvar = dvar * (tw[1] - tw[0])
    elif dv == "Amp":
        dvar = load_mat_variable(find_by_suffix("AmpMatx.mat"), "AmpMatx")
    elif dv == "Vel":
        dvar = load_mat_variable(find_by_suffix("VelMatx.mat"), "VelMatx")
    elif dv in {"Ang", "Cos", "Sin"}:
        dvar = load_mat_variable(find_by_suffix("AngMatx.mat"), "AngMatx")
    else:
        # Fallback to beh/info fields
        beh = load_mat_variable(find_by_suffix("beh.mat"), "behData")
        info = load_mat_variable(find_by_suffix("info.mat"), "subInfo")
        if hasattr(beh, "dtype") and beh.dtype.names and dv in beh.dtype.names:
            dvar = beh[dv].flatten()
        elif (hasattr(info, dv)) or (isinstance(info, dict) and dv in info):
            dvar = getattr(info, dv) if hasattr(info, dv) else info[dv]
        else:
            raise ValueError(f"Unsupported DV: {dv}")

    # Time-window averaging may require info.epoch
    if tw is not None:
        try:
            info = load_mat_variable(find_by_suffix("info.mat"), "subInfo")
            epoch = getattr(info, "epoch") if hasattr(info, "epoch") else (info.get("epoch") if isinstance(info, dict) else None)
            if epoch is not None:
                tw_shifted = [tw[0] - epoch[0], tw[1] - epoch[0]]
                dvar = np.nanmean(dvar[tw_shifted[0]:tw_shifted[1], :], axis=0)
        except Exception:
            pass

    dvar = np.atleast_2d(dvar)
    if return_meta:
        # If explicit meta not loaded, best-effort to provide it
        if beh_meta is None:
            try:
                beh_meta = load_mat_variable(find_by_suffix("beh.mat"), "behData")
            except Exception:
                beh_meta = None
        if info_meta is None:
            try:
                info_meta = load_mat_variable(find_by_suffix("info.mat"), "subInfo")
            except Exception:
                info_meta = None
        return dvar, beh_meta, info_meta
    return dvar


def get_iv(subj_dir, iv, iv_val):
    """
    Extract independent variable masks and labels.
    """
    # find files by suffix to match naming pattern
    def find_by_suffix(suffix):
        files = [f for f in os.scandir(subj_dir) if f.is_file() and f.name.endswith(suffix)]
        if not files:
            raise FileNotFoundError(f"No file ending with {suffix} found in {subj_dir}")
        return files[0].path

    def load_mat_variable(path, var):
        try:
            d = loadmat(path, squeeze_me=True, struct_as_record=False)
            data_keys = [k for k in d.keys() if not k.startswith("__")]
            if var in d:
                return d[var]
            if data_keys:
                return d[data_keys[0]]
        except NotImplementedError:
            with h5py.File(path, "r") as f:
                if var in f:
                    obj = f[var]
                    if isinstance(obj, h5py.Dataset):
                        return obj[()]
                    if isinstance(obj, h5py.Group):
                        return {k: (v[()] if isinstance(v, h5py.Dataset) else v) for k, v in obj.items()}
                for k, v in f.items():
                    if isinstance(v, h5py.Dataset):
                        return v[()]
                    if isinstance(v, h5py.Group):
                        return {kk: (vv[()] if isinstance(vv, h5py.Dataset) else vv) for kk, vv in v.items()}
        raise KeyError(f"Variable {var} not found in {path}")

    beh = load_mat_variable(find_by_suffix("beh.mat"), "behData")
    info = load_mat_variable(find_by_suffix("info.mat"), "subInfo")
    rmflag = load_mat_variable(find_by_suffix("rmFlag.mat"), "rmFlag")

    def get_field(container, name):
        # Support dict, scipy mat_struct, numpy structured arrays
        if isinstance(container, dict):
            return container[name]
        if hasattr(container, name):
            return getattr(container, name)
        if hasattr(container, "dtype") and container.dtype.names and name in container.dtype.names:
            return container[name]
        raise TypeError("Unsupported container for field access")

    # determine number of trials
    if hasattr(beh, "shape"):
        num_trial = beh.shape[0]
    elif isinstance(beh, dict):
        first_key = next(iter(beh))
        num_trial = np.asarray(beh[first_key]).shape[0]
    else:
        raise ValueError("Unsupported beh format")

    if iv is None or iv == "none":
        overall = np.asarray(get_field(rmflag, "overall")).flatten()
        rm_idx = [~overall]
        labels = pd.DataFrame()
    else:
        if iv_val is None:
            alias_map = {"group": "group_id", "cond": "condition", "conditions": "condition"}

            def resolve_name_in_info(info_obj, name):
                # try exact
                if hasattr(info_obj, name):
                    return name
                # alias
                if name in alias_map and hasattr(info_obj, alias_map[name]):
                    return alias_map[name]
                # case-insensitive match over attrs
                attrs = [k for k in dir(info_obj) if not k.startswith('_')]
                for a in attrs:
                    if a.lower() == name.lower():
                        return a
                return None

            def resolve_name_in_beh(beh_obj, name):
                # exact
                if hasattr(beh_obj, "dtype") and beh_obj.dtype.names and name in beh_obj.dtype.names:
                    return name
                if isinstance(beh_obj, dict) and name in beh_obj:
                    return name
                # alias
                if name in alias_map:
                    alias = alias_map[name]
                    if hasattr(beh_obj, "dtype") and beh_obj.dtype.names and alias in beh_obj.dtype.names:
                        return alias
                    if isinstance(beh_obj, dict) and alias in beh_obj:
                        return alias
                # case-insensitive
                names = []
                if hasattr(beh_obj, "dtype") and beh_obj.dtype.names:
                    names = list(beh_obj.dtype.names)
                elif isinstance(beh_obj, dict):
                    names = list(beh_obj.keys())
                for n in names:
                    if n.lower() == name.lower():
                        return n
                # substring match
                for n in names:
                    if name.lower() in n.lower():
                        return n
                return None

            def get_series(container, col):
                # Try from beh first, then from info
                bname = resolve_name_in_beh(container, col)
                if bname is not None:
                    if hasattr(container, "dtype") and container.dtype.names:
                        return container[bname].flatten()
                    return np.asarray(container[bname]).flatten()
                return None

            data = {}
            for col in iv:
                s = get_series(beh, col)
                if s is None:
                    # try from info
                    if isinstance(info, dict):
                        # try exact, alias, ci, substring
                        key = col
                        if key not in info and col in alias_map and alias_map[col] in info:
                            key = alias_map[col]
                        if key not in info:
                            # case-insensitive
                            ci = [k for k in info.keys() if k.lower() == col.lower()]
                            if ci:
                                key = ci[0]
                        if key not in info:
                            # substring
                            sub = [k for k in info.keys() if col.lower() in k.lower()]
                            if sub:
                                key = sub[0]
                        if key in info:
                            s = np.asarray(info[key]).flatten()
                        else:
                            s = None
                    else:
                        actual = resolve_name_in_info(info, col)
                        if actual is not None:
                            s = np.asarray(getattr(info, actual)).flatten()
                        else:
                            raise ValueError(f"IV column '{col}' not found in beh or info")
                data[col] = s
            groupings = pd.DataFrame(data)
            group_idx = pd.Series(list(map(tuple, groupings.values))).astype("category").cat.codes
            n_cond = group_idx.max() + 1
            rm_idx = []
            labels = []
            for i in range(n_cond):
                mask = (group_idx == i).values
                overall = np.asarray(get_field(rmflag, "overall")).flatten()
                rm_idx.append((~overall) & mask)
                labels.append(groupings.iloc[i:i+1])
            labels = pd.concat(labels, ignore_index=True)
        else:
            rm_idx = []
            for row in iv_val:
                mask = np.ones(num_trial, dtype=bool)
                for col, val in zip(iv, row):
                    bname = resolve_name_in_beh(beh, col)
                    if bname is not None:
                        series = beh[bname].flatten() if hasattr(beh, "dtype") and beh.dtype.names else np.asarray(beh[bname]).flatten()
                    elif isinstance(info, dict) and col in info:
                        series = np.asarray(info[col]).flatten()
                    else:
                        actual = resolve_name_in_info(info, col)
                        if actual is not None:
                            series = np.asarray(getattr(info, actual)).flatten()
                        else:
                            raise ValueError(f"IV column '{col}' not found in beh or info")
                    mask &= (series == val)
                overall = np.asarray(get_field(rmflag, "overall")).flatten()
                rm_idx.append((~overall) & mask)
            labels = pd.DataFrame(iv_val, columns=iv)

    return labels, rm_idx
