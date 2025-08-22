# Proposal: Return Z Beta Estimates in LSS Functions

## Current Situation

After examining all LSS implementations, here's what happens to the Z (experimental regressor) beta estimates:

### What All Implementations Currently Do:
1. **Combine Z with nuisance regressors** into a single confound matrix: `X_confounds = cbind(Z, Nuisance)`
2. **Project out the combined confounds** from both Y and trial regressors (X)
3. **Compute trial-wise betas** on the projected data
4. **Return only trial betas** (T × V matrix)
5. **Discard the Z beta estimates** that were computed during the projection step

### Key Finding:
The Z beta estimates ARE computed during the projection step in all implementations, but they're thrown away. We need to capture and return them.

## Proposed Solution

### 1. Modify the Projection Step
Instead of treating Z and Nuisance as a single confound matrix, we need to:

```r
# Current approach:
X_confounds <- cbind(Z, Nuisance)  # Combine everything
# Project out all confounds together

# Proposed approach:
# Step 1: Project out Nuisance from Y, X, and Z
if (!is.null(Nuisance)) {
  Y_clean <- residualize(Y, Nuisance)
  X_clean <- residualize(X, Nuisance) 
  Z_clean <- residualize(Z, Nuisance)
} else {
  Y_clean <- Y
  X_clean <- X
  Z_clean <- Z
}

# Step 2: Fit Z to clean data and get Z betas
Z_betas <- solve(t(Z_clean) %*% Z_clean) %*% t(Z_clean) %*% Y_clean

# Step 3: Project out Z from Y and X for trial estimation
Y_final <- Y_clean - Z_clean %*% Z_betas
X_final <- residualize(X_clean, Z_clean)

# Step 4: Compute trial betas on final projected data
trial_betas <- lss_computation(Y_final, X_final)
```

### 2. Return Structure
Change the return value from a simple matrix to a list:

```r
# Current return:
matrix(trial_betas)  # T × V

# Proposed return:
list(
  trial_betas = matrix(trial_betas),    # T × V (trial estimates)
  Z_betas = matrix(Z_betas),            # F × V (experimental regressor estimates)
  trials = rownames(trial_betas),       # Trial names
  voxels = colnames(trial_betas),       # Voxel names
  regressors = colnames(Z)              # Z regressor names
)
```

### 3. Backward Compatibility
Provide a parameter to control return format:

```r
lss(..., return_format = c("matrix", "list"))

# return_format = "matrix": current behavior (trial betas only)
# return_format = "list": new behavior (both trial and Z betas)
```

## Implementation Plan

### Phase 1: Core Infrastructure
1. **Create helper functions** for the new projection workflow
2. **Update the main lss() function** to support the new return format
3. **Modify .lss_r_optimized()** as the reference implementation

### Phase 2: Update All Implementations
1. **Update .lss_naive()** - easiest to modify and test
2. **Update .lss_r_vectorized()** and **lss_fast()**
3. **Update .lss_cpp()** and related C++ functions
4. **Update .lss_cpp_optimized()** and **lss_fused_optim_cpp()**

### Phase 3: C++ Implementation Updates
1. **Modify compute_residuals_cpp()** to handle separate projection steps
2. **Create new C++ functions** for Z beta computation
3. **Update lss_fused_optim_cpp()** for the new workflow

### Phase 4: Testing and Documentation
1. **Create comprehensive tests** for Z beta estimates
2. **Update documentation** and examples
3. **Add vignette section** on interpreting Z betas

## Specific Implementation Details

### Helper Functions Needed:

```r
# Residualize Y against X: Y - X %*% solve(t(X) %*% X) %*% t(X) %*% Y
.residualize <- function(Y, X) {
  if (is.null(X)) return(Y)
  XtX_inv <- chol2inv(chol(crossprod(X)))
  Y - X %*% (XtX_inv %*% crossprod(X, Y))
}

# Compute beta estimates: solve(t(X) %*% X) %*% t(X) %*% Y
.compute_betas <- function(Y, X) {
  XtX_inv <- chol2inv(chol(crossprod(X)))
  XtX_inv %*% crossprod(X, Y)
}

# Format return value based on return_format
.format_lss_output <- function(trial_betas, Z_betas = NULL, 
                               return_format = "matrix", ...) {
  if (return_format == "matrix") {
    return(trial_betas)
  } else {
    return(list(
      trial_betas = trial_betas,
      Z_betas = Z_betas,
      ...
    ))
  }
}
```

### C++ Functions Needed:

```cpp
// Separate projection for nuisance vs experimental regressors
List project_separate_cpp(const arma::mat& Y,
                          const arma::mat& X_trials,
                          const arma::mat& Z_experimental,
                          const arma::mat& Nuisance);

// Compute beta estimates
arma::mat compute_betas_cpp(const arma::mat& Y, const arma::mat& X);
```

## Benefits

1. **Complete Information**: Users get all beta estimates computed during LSS
2. **Experimental Design Insights**: Can examine condition effects, block effects, etc.
3. **Debugging and Validation**: Can verify that nuisance regression worked properly
4. **Flexibility**: Users can choose what information they need
5. **Backward Compatibility**: Existing code continues to work unchanged

## Considerations

1. **Memory Usage**: Returning Z betas increases memory usage, but typically Z has few columns
2. **API Consistency**: Need to decide on the best return structure
3. **Performance**: The new projection workflow may be slightly slower but more informative
4. **Documentation**: Need clear examples of when and how to use Z betas

## Recommended First Steps

1. **Implement the helper functions** 
2. **Modify .lss_r_optimized() as a prototype** with `return_format` parameter
3. **Create comprehensive tests** comparing old vs new approaches
4. **Get user feedback** on the proposed return structure
5. **Roll out to other implementations** once the design is validated 