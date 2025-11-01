# fmrilss Command-Line Scripts

This directory contains command-line tools for trial-wise beta estimation using fmrilss.

## Scripts Overview

### estimate_lss_bids.R ⭐ (Recommended)

**Full BIDS-integrated pipeline** combining automatic dataset navigation with modern fmrilss estimation.

**Key features:**
- ✅ Automatic BIDS dataset discovery via `bidser`
- ✅ Confound extraction with PCA variance retention
- ✅ Multi-run and multi-session support
- ✅ Multiple estimation backends (r, cpp_optimized, oasis)
- ✅ Optional AR prewhitening
- ✅ Polynomial detrending (--polort)
- ✅ No AFNI dependency

**Basic usage:**
```bash
./estimate_lss_bids.R \
  --bids_path=/data/myproject \
  --subid=01 \
  --task=stroop \
  --tr=2.0 \
  --method=oasis \
  --outdir=betas
```

**With confounds:**
```bash
./estimate_lss_bids.R \
  --bids_path=/data/myproject \
  --subid=01 \
  --task=stroop \
  --tr=2.0 \
  --confounds=confounds.txt \
  --percvar=95 \
  --method=cpp_optimized \
  --prewhiten=ar1
```

**With mask and session:**
```bash
./estimate_lss_bids.R \
  --bids_path=/data/myproject \
  --subid=01 \
  --task=faces \
  --bids_session=01 \
  --tr=2.0 \
  --mask=task-faces_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz \
  --method=oasis
```

**Concatenate runs:**
```bash
./estimate_lss_bids.R \
  --bids_path=/data/myproject \
  --subid=01 \
  --task=nback \
  --tr=2.0 \
  --concatenate=true \
  --out=nback_betas \
  --method=oasis
```

### estimate_lss.R

**Standalone script** for single-file processing (no BIDS integration).

**Best for:**
- Processing individual scan files
- Non-BIDS datasets
- Quick testing

**Usage:**
```bash
./estimate_lss.R \
  --bold=/path/to/bold.nii.gz \
  --events=/path/to/events.tsv \
  --tr=2.0 \
  --method=oasis \
  --outdir=results
```

### estimate_betas.R

**Legacy AFNI-based script** (requires AFNI installation).

**Note:** Use `estimate_lss_bids.R` instead for modern workflows without AFNI dependency.

## Parameters Reference

### BIDS Parameters (estimate_lss_bids.R only)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bids_path` | `.` | Path to BIDS dataset root |
| `--subid` | *required* | Subject ID (e.g., `01`) |
| `--task` | *required* | Task name (e.g., `stroop`) |
| `--bids_session` | `""` | Session ID (optional, e.g., `01`) |
| `--deriv_folder` | `derivatives/fmriprep` | Preprocessing derivatives folder |
| `--space` | `MNI152NLin2009cAsym` | BOLD image space |

### Design Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tr` | *required* | Repetition time in seconds |
| `--basis` | `spmg1` | HRF basis (`spmg1` or `spmg3`) |
| `--duration` | `0` | Default event duration (seconds) |
| `--onset_col` | `onset` | Column name for onsets in events.tsv |
| `--duration_col` | `NULL` | Column name for durations (optional) |
| `--milliseconds` | `false` | Whether onsets are in milliseconds |
| `--polort` | `3` | Polynomial detrending order |

### Confound Parameters (estimate_lss_bids.R only)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--confounds` | `NULL` | Text file listing confound variable names |
| `--percvar` | `95` | Percentage of confound variance to retain via PCA |

**Example confounds.txt:**
```
trans_x
trans_y
trans_z
rot_x
rot_y
rot_z
framewise_displacement
```

### Estimation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--method` | `oasis` | Estimation method: `r`, `r_optimized`, `cpp_optimized`, `oasis`, `naive` |
| `--prewhiten` | `none` | Prewhitening: `none`, `ar1`, `auto` |
| `--mask` | `NULL` | Binary mask file (auto-mask if not provided) |

### Output Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `betas` | Output directory name |
| `--out` | `betas` | Output file stem |
| `--concatenate` | `false` | Concatenate all runs into single file |
| `--write-nifti` | `true` | Write 4D NIfTI output |
| `--write-csv` | `false` | Write CSV output (trials × voxels) |

## Estimation Methods

### r_optimized (Default for non-OASIS)
- Pure R implementation with optimizations
- Good balance of speed and compatibility
- No compilation required

### cpp_optimized
- Parallelized C++ with OpenMP
- **Fastest** for large datasets
- Requires Rcpp compilation

### oasis
- Optimized Analytic Single-pass Inverse Solution
- Supports ridge regularization
- Mathematically equivalent to LSS
- Good for high-collinearity designs
- **Recommended for most use cases**

### naive
- Simple reference implementation
- Use only for testing

## Prewhitening Options

### none (Default)
- No temporal autocorrelation correction
- Fastest but assumes independent errors

### ar1
- First-order autoregressive correction
- Good for typical fMRI noise

### auto
- Automatic AR order selection (up to AR(4))
- Most robust but slower
- Uses AIC/BIC model selection

## Output Files

### NIfTI Format (4D)
- Dimensions: X × Y × Z × ntrials
- Each volume = beta map for one trial
- Preserves spatial structure

### CSV Format (optional)
- Dimensions: ntrials × voxels
- Columns: voxel indices
- Rows: trial betas
- Good for quick inspection

## Examples

### Example 1: Basic BIDS Processing
```bash
./estimate_lss_bids.R \
  --bids_path=/data/studyXYZ \
  --subid=03 \
  --task=workingmemory \
  --tr=2.5 \
  --method=oasis \
  --outdir=lss_betas
```

### Example 2: With Motion Confounds
Create `motion_confounds.txt`:
```
trans_x
trans_y
trans_z
rot_x
rot_y
rot_z
framewise_displacement
```

Run:
```bash
./estimate_lss_bids.R \
  --bids_path=/data/studyXYZ \
  --subid=03 \
  --task=workingmemory \
  --tr=2.5 \
  --confounds=motion_confounds.txt \
  --percvar=95 \
  --method=cpp_optimized \
  --prewhiten=ar1
```

### Example 3: Multi-Session with Custom Mask
```bash
./estimate_lss_bids.R \
  --bids_path=/data/studyXYZ \
  --subid=03 \
  --task=faces \
  --bids_session=wave2 \
  --tr=2.0 \
  --mask=task-faces_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz \
  --method=oasis \
  --concatenate=true
```

### Example 4: Advanced Options
```bash
./estimate_lss_bids.R \
  --bids_path=/data/studyXYZ \
  --subid=03 \
  --task=gambling \
  --tr=2.0 \
  --basis=spmg3 \
  --duration=1.5 \
  --confounds=full_confounds.txt \
  --percvar=98 \
  --method=oasis \
  --prewhiten=auto \
  --polort=4 \
  --concatenate=true \
  --write-csv=true \
  --outdir=gambling_lss
```

## Troubleshooting

### "bidser package required"
Install bidser from GitHub:
```r
remotes::install_github("bbuchsbaum/bidser")
```

### "No preprocessed scans found"
Check:
1. BIDS path is correct
2. fMRIPrep derivatives exist in `derivatives/fmriprep`
3. Subject ID matches (without `sub-` prefix)
4. Task name matches (without `task-` prefix)

### "Design rows != BOLD timepoints"
Event onsets may be incorrect:
- Check if `--milliseconds=true` is needed
- Verify onset column name with `--onset_col`
- Ensure TR is correct

### "Mask not found"
Provide full mask filename including BIDS entities:
```bash
--mask=sub-01_task-faces_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz
```

### Slow estimation
Try:
1. Use `--method=cpp_optimized` for speed
2. Reduce confound PCA: `--percvar=90`
3. Skip prewhitening: `--prewhiten=none`
4. Use a more restrictive mask

## Performance Tips

1. **Use cpp_optimized** for datasets with >10,000 voxels
2. **OASIS method** is memory-efficient for long time series
3. **Concatenate runs** (`--concatenate=true`) for across-run analyses
4. **PCA confounds** at 90-95% captures most variance efficiently
5. **AR(1) prewhitening** is usually sufficient; avoid `auto` unless needed

## Dependencies

### Required R packages:
- fmrilss
- fmrihrf
- bidser (BIDS integration only)
- RNifti

### Optional:
- fmridesign (better HRF handling)
- fmriAR (prewhitening)
- data.table (faster I/O)

## Citation

If you use fmrilss in your research, please cite:

> Buchsbaum, B.R. (2024). fmrilss: Efficient Least Squares Separate analysis for fMRI. R package version X.X.X. https://github.com/bbuchsbaum/fmrilss
