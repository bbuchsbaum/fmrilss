

### **Granular Proposal: Voxel-wise HRF Estimation in `fmrilss`**

**1. Executive Summary**

This document specifies the design and implementation of an extension to the `fmrilss` package to support voxel-wise Hemodynamic Response Function (HRF) estimation within the Least Squares-Separate (LSS) framework. The design prioritizes computational efficiency, API clarity, and seamless integration with the existing `fmrilss` and `fmrihrf` packages.

**2. Core Philosophy & Design Principles**

The implementation will adhere to three core principles:

1.  **Estimator, Not Manager:** The new functions will operate on numeric matrices (`Y`, `Z`) and event tables, maintaining `fmrilss`'s role as a statistical estimator, not an image data manager.
2.  **Modularity:** The workflow is split into two distinct, user-inspectable stages: (1) HRF coefficient estimation and (2) trial-wise beta estimation using those HRFs. This enhances debuggability and aligns with established neuroimaging analysis patterns.
3.  **Tool Re-use:** The `fmrihrf` package will be leveraged for all HRF basis generation, convolution, and reconstruction tasks, ensuring consistency and minimizing code duplication.

**3. API Specification**

We will introduce two new exported functions.

#### **3.1 `estimate_voxel_hrf()` - Stage 1: HRF Estimation**

This function performs a multi-response GLM to estimate the HRF basis coefficients for every voxel simultaneously.

```R
#' @title Estimate Voxel-wise HRF Basis Coefficients
#' @description Fits a GLM to estimate the coefficients for a given HRF basis set for each voxel.
#'
#' @param Y A numeric matrix of BOLD data (time x voxels).
#' @param events A data.frame with `onset`, `duration`, and `condition` columns.
#' @param basis An `HRF` object from the `fmrihrf` package (e.g., `HRF_FIR`, `HRF_BSPLINE`).
#'   This defines the basis set used for HRF estimation.
#' @param nuisance_regs A numeric matrix of nuisance regressors (time x num_nuisance).
#'
#' @return An S3 object of class `VoxelHRF`, which is a list containing:
#'   \item{coefficients}{A matrix of HRF basis coefficients (num_params x num_voxels).}
#'   \item{basis}{The `HRF` object used for estimation.}
#'   \item{conditions}{A character vector of the conditions modeled.}
#'   \item{...}{Other metadata for provenance.}
#' @export
estimate_voxel_hrf <- function(Y, events, basis, nuisance_regs = NULL) {
  # Implementation...
}
```

#### **3.2 `lss_with_hrf()` - Stage 2: LSS Beta Estimation**

This function takes the estimated HRF coefficients and performs the LSS analysis, generating the trial design matrix on-the-fly for each voxel.

```R
#' @title Perform LSS using Voxel-wise HRFs
#' @description Computes trial-wise beta estimates using voxel-specific HRFs.
#'
#' @param Y A numeric matrix of BOLD data (time x voxels).
#' @param events A data.frame with `onset`, `duration`, and `condition` columns.
#' @param hrf_estimates A `VoxelHRF` object returned by `estimate_voxel_hrf()`.
#' @param nuisance_regs A numeric matrix of nuisance regressors.
#' @param engine A character string specifying the computational engine: "C++" (default) or "R".
#' @param chunk_size An integer specifying the number of voxels to process per batch (C++ engine only).
#' @param verbose A logical indicating whether to display a progress bar.
#'
#' @return An S3 object of class `LSSBeta`, which is a list containing:
#'   \item{betas}{A matrix of LSS beta estimates (num_trials x num_voxels), potentially file-backed.}
#'   \item{...}{Other metadata.}
#' @export
lss_with_hrf <- function(Y, events, hrf_estimates, nuisance_regs = NULL,
                         engine = "C++", chunk_size = 5000, verbose = TRUE) {
  # Implementation...
}
```

**4. Computational Workflow**

#### **4.1 `estimate_voxel_hrf()` Workflow**

This stage is a straightforward multi-response GLM.

1.  **Design Matrix Construction (R):**
    *   Use `fmrihrf::regressor_set()` with the provided `basis` object to generate the full rank design matrix `X_basis` for all conditions and basis functions.
    *   `cbind` `X_basis` with `nuisance_regs` to form the final design matrix `X_full`.
2.  **GLM Solve (C++/Armadillo):**
    *   The core is a single, highly efficient call: `B = arma::solve(X_full, Y)`.
    *   This solves for all voxel coefficients `B` in one operation.
3.  **Return Structure (R):**
    *   Package the relevant coefficient rows from `B` and metadata into the `VoxelHRF` object.

#### **4.2 `lss_with_hrf()` Workflow**

This is the core computational engine. The R function will act as a wrapper for the C++ engine.

1.  **R-side Preparation:**
    *   Validate inputs.
    *   Extract event timing (`onset`, `duration`) and condition labels into simple vectors.
    *   Evaluate the HRF basis functions on a fine time grid to get `hrf_basis_kernels` (matrix: `fine_time x num_basis`). This is cached.
    *   Cache event onset indices relative to the fine time grid.
    *   If `bigmemory` is used, initialize the file-backed output matrix for `betas`.
    *   Initialize progress bar if `verbose=TRUE`.
    *   Call the C++ engine, passing pointers/references to all data.

2.  **C++ Engine (`lss_engine_vox_hrf`):**
    *   The engine will loop through `Y` and `hrf_estimates$coefficients` in chunks.
    *   An outer `OpenMP` parallel pragma will be placed on the loop over voxels within a chunk.
    *   **Inside the parallel loop (per-voxel `v`):**
        1.  **Thread-Private Storage:** Initialize a private vector to store the beta estimates for the current voxel.
        2.  **On-the-fly HRF Kernel Reconstruction:** Compute the voxel-specific HRF kernel: `arma::vec hrf_v_kernel = hrf_basis_kernels * hrf_coefficients.col(v);`
        3.  **On-the-fly Trial Regressor Generation:** Using the cached event indices, convolve `hrf_v_kernel` with the event timing to produce the trial-wise design matrix `C_v` (`arma::mat`). This will use `arma::conv`.
        4.  **LSS Solve:** Feed `Y.col(v)`, `C_v`, and `nuisance_regs` into the existing, stable LSS solver logic from `fmrilss`.
        5.  **Numerical Stability:** The LSS solver's `XtX` calculation will be regularized with a small ridge term (`1e-8 * arma::eye`) before `arma::solve` to ensure stability.
        6.  Store the resulting betas in the thread-private vector.
    *   **After the parallel loop (per chunk):**
        *   An `#pragma omp critical` block will be used to safely write the results from each thread's private storage into the shared output matrix (e.g., the memory-mapped `bigmemory` file).
    *   Update the progress bar from the R side after each chunk is processed.

**5. Data Handling & Memory Management**

To handle whole-brain datasets, the `betas` and `coefficients` matrices will be backed by the `bigmemory` package. The C++ code will receive the file path or pointer to the memory-mapped object, allowing it to write results out-of-core and avoid exhausting RAM. The `DelayedArray` suggestion is noted as a valuable future enhancement for Bioconductor interoperability.

**6. Validation & Testing Strategy**

The test suite will be comprehensive and modular:

1.  **Unit Tests (`testthat`):**
    *   **`estimate_voxel_hrf()`:**
        *   Verify correct design matrix construction.
        *   *Synthetic Recovery:* Generate data with known HRF coefficients and confirm their recovery.
    *   **`lss_with_hrf()`:**
        *   *Synthetic Recovery:* Generate data with known voxel-wise HRFs and trial betas; confirm recovery of the betas.
        *   *Equivalence Test:* When all voxels are given identical (e.g., canonical) HRF coefficients, results *must* match the output of the standard `fmrilss::lss()` function.
        *   *Basis Aliasing Test:* Generate data with a known spline HRF, fit with a more flexible FIR basis, then use `hrf_from_coefficients` to reconstruct the shape. The reconstructed shape must match the ground-truth spline.
2.  **Integration Test:** A full-pipeline test using a small, simulated dataset.
3.  **Real Data Validation:** Replicate key findings from published literature (e.g., latency differences in HCP retinotopy data) to ensure real-world validity.

**7. Documentation & User Experience**

*   **Vignette:** A new vignette titled "Advanced Modeling: Voxel-wise HRFs" will be created. It will demonstrate the full workflow: simulating data, estimating HRFs, performing visual QC on the HRF maps, running LSS, and interpreting the results.
*   **Function Documentation:** The R-level functions will be thoroughly documented with examples.
*   **Progress Reporting:** The `verbose = TRUE` argument in `lss_with_hrf` will be implemented using the `progress` package or similar, providing crucial feedback for long-running analyses.

**8. Next Steps**

1.  **Branch:** Create a `feature/voxel-hrf` branch.
2.  **Stubbing:** Implement the R function stubs for `estimate_voxel_hrf` and `lss_with_hrf` with input validation.
3.  **C++ Implementation:**
    *   Develop the C++ kernel for `estimate_voxel_hrf`.
    *   Develop the core `lss_engine_vox_hrf` C++ kernel for `lss_with_hrf`, focusing on the parallelized on-the-fly convolution and LSS solve.
4.  **Unit Testing:** Write `testthat` scripts for each function as it is developed, starting with synthetic data recovery.
5.  **Benchmarking:** Profile the C++ engine on a medium-sized dataset (~100k voxels, 600 TRs) to ensure performance goals are met and identify any bottlenecks.

This plan provides a clear, robust, and efficient path to implementing a powerful and much-needed feature in the `fmrilss` package.

### **Sprint 1: The Foundation - HRF Estimation & API Scaffolding**

**Goal:** Implement the complete `estimate_voxel_hrf` workflow and set up the API and data preparation logic for `lss_with_hrf`. By the end of this sprint, a user can estimate voxel-wise HRFs and have a clear (though not yet functional) interface for the next step.

---

**Ticket S1-T1: Project Setup & API Definition**

*   **Description:** Initialize the development environment. Create a new feature branch `feature/voxel-hrf`. Create placeholder files for new R functions (`R/voxel_hrf.R`) and C++ code (`src/voxel_hrf.cpp`). Define the S3 class structures for `VoxelHRF` and `LSSBeta` objects.
*   **AC:**
    1.  The `feature/voxel-hrf` branch exists in the repository.
    2.  `R/voxel_hrf.R` is created and contains the R function stubs for `estimate_voxel_hrf()` and `lss_with_hrf()`.
    3.  `src/voxel_hrf.cpp` is created.
    4.  Initial S3 class definitions are documented (e.g., using roxygen2 `@return` tags).

**Ticket S1-T2: Implement R-side Logic for `estimate_voxel_hrf`**

*   **Description:** Implement the R portion of `estimate_voxel_hrf`. This involves validating all inputs (`Y`, `events`, `basis`, `nuisance_regs`) and constructing the full design matrix `X_full` by leveraging `fmrihrf::regressor_set` for the basis functions and `cbind` for nuisance regressors. The function should then call a placeholder C++ function.
*   **AC:**
    1.  The function robustly validates all input types and dimensions.
    2.  It correctly uses `fmrihrf` to generate the `X_basis` matrix.
    3.  It correctly combines `X_basis` and `nuisance_regs`.
    4.  The prepared `X_full` and `Y` matrices are passed to a C++ function stub.

**Ticket S1-T3: Develop C++ Kernel for HRF Estimation**

*   **Description:** Implement the C++ backend for `estimate_voxel_hrf`. This is a single, efficient GLM solver. The function will accept the `X_full` and `Y` matrices (as `arma::mat`) and return the coefficient matrix `B` using `arma::solve`.
*   **AC:**
    1.  A C++ function `estimate_hrf_cpp(const arma::mat& X, const arma::mat& Y)` is created.
    2.  The function returns `arma::solve(X, Y)`.
    3.  The function is exposed to R via `Rcpp::export`.

**Ticket S1-T4: Integrate and Finalize `estimate_voxel_hrf`**

*   **Description:** Connect the R function from S1-T2 to the C++ kernel from S1-T3. The R function will receive the full coefficient matrix `B` from C++, extract the relevant rows corresponding to the HRF basis functions, and package the results into a final `VoxelHRF` S3 object with all specified metadata.
*   **AC:**
    1.  `estimate_voxel_hrf()` runs end-to-end without errors on a test dataset.
    2.  The returned object is of class `VoxelHRF` and contains `coefficients`, `basis`, and `conditions`.
    3.  The dimensions of the `coefficients` matrix are correct (`num_params x num_voxels`).

**Ticket S1-T5: Unit Tests for `estimate_voxel_hrf`**

*   **Description:** Create a new test file `tests/testthat/test-voxel-hrf.R`. Write comprehensive unit tests for `estimate_voxel_hrf`. This includes tests for input validation and, most importantly, a synthetic data recovery test.
*   **AC:**
    1.  A test that generates `Y` from known HRF coefficients successfully recovers those coefficients within a small tolerance (e.g., RMSE < 1e-6).
    2.  Tests for invalid inputs (e.g., mismatched dimensions, incorrect object types) correctly throw errors.
    3.  All tests pass in `devtools::test()`.

**Ticket S1-T6: Implement R-side Data Preparation for `lss_with_hrf`**

*   **Description:** In the `lss_with_hrf` R function, implement the data preparation logic that occurs *before* calling the C++ engine. This includes validating inputs, extracting event timings into vectors, evaluating the HRF basis kernels on a fine time grid (`hrf_basis_kernels`), and caching event onset indices.
*   **AC:**
    1.  The function validates the `hrf_estimates` object structure.
    2.  It correctly generates the `hrf_basis_kernels` matrix.
    3.  It correctly prepares event data (onsets, durations) as simple vectors.
    4.  The function currently ends by calling a placeholder C++ function stub with the prepared data structures.

---

### **Sprint 2: Core LSS Engine, Integration, and Polishing**

**Goal:** Implement the high-performance C++ engine for `lss_with_hrf`, integrate it with the R wrapper, add performance and usability features, and complete all documentation. By the end of this sprint, the entire feature is production-ready.

---

**Ticket S2-T1: Develop C++ LSS Engine - Serial Voxel Loop**

*   **Description:** Implement the core logic for the `lss_with_hrf` C++ engine (`lss_engine_vox_hrf`). Focus on correctness first with a **serial** loop over voxels. For each voxel: reconstruct its HRF kernel (`hrf_basis_kernels * coefficients_v`), generate the trial design matrix `C_v` via on-the-fly convolution (`arma::conv`), and perform the LSS solve using existing `fmrilss` logic.
*   **AC:**
    1.  The C++ function can process a single voxel and produce beta estimates.
    2.  For a single voxel, the output exactly matches a pure R implementation of the same logic.
    3.  The LSS solve includes the `1e-8 * I` ridge regularization for stability.

**Ticket S2-T2: Parallelize C++ Engine with OpenMP**

*   **Description:** Introduce `OpenMP` to parallelize the voxel loop from S2-T1. Implement chunking logic. Ensure thread safety by having each thread write its results to a private vector/matrix and then copy them to the shared output matrix within a `#pragma omp critical` block after the parallel region for that chunk.
*   **AC:**
    1.  The parallelized engine produces bit-for-bit identical results to the serial engine from S2-T1.
    2.  Performance scales with the number of available cores (verified via benchmarking on a multi-core machine).
    3.  The engine correctly processes data in chunks.

**Ticket S2-T3: Integrate `bigmemory` and Progress Reporting**

*   **Description:** Modify the C++ engine to write its output directly to a memory-mapped file provided via `bigmemory` (passed from R as an `XPtr`). Implement a callback mechanism or simple counter so the R wrapper can query progress and update a `progress` bar after each chunk.
*   **AC:**
    1.  The C++ function can successfully receive a `bigmemory` pointer and write results to the file-backed matrix.
    2.  The R function `lss_with_hrf` can display a functional progress bar when `verbose=TRUE`.
    3.  The final `LSSBeta` object correctly references the file-backed `betas` matrix.

**Ticket S2-T4: Unit Tests for `lss_with_hrf`**

*   **Description:** Add comprehensive unit tests for `lss_with_hrf`. This must include a synthetic data recovery test and the critical equivalence test.
*   **AC:**
    1.  A test that generates data with known voxel-wise HRFs and trial betas successfully recovers the betas.
    2.  **Equivalence Test:** A test where all voxels are assigned the *same* canonical HRF coefficients must produce results that are numerically identical to the existing `fmrilss::lss()` function.
    3.  **Basis Aliasing Test:** A test that fits a flexible FIR basis to data generated from a known spline HRF, then reconstructs the HRF shape, confirms the shape matches the ground truth.
    4.  All tests pass.

**Ticket S2-T5: Create Vignette**

*   **Description:** Write a new package vignette (`voxel-wise-hrf.Rmd`) that walks a user through the entire two-step workflow. The vignette should include code and visualizations for: (1) Simulating data, (2) Running `estimate_voxel_hrf`, (3) Plotting example HRF shapes and a map of a parameter (e.g., time-to-peak), (4) Running `lss_with_hrf`, and (5) Analyzing the final trial betas.
*   **AC:**
    1.  The vignette is written and renders correctly to HTML.
    2.  All code chunks are executable and produce the documented output/plots.
    3.  The narrative clearly explains the motivation and steps.

**Ticket S2-T6: Finalize All Documentation**

*   **Description:** Write complete Roxygen2 documentation for `estimate_voxel_hrf`, `lss_with_hrf`, and the new S3 objects (`VoxelHRF`, `LSSBeta`). Ensure all parameters are described and each function has a runnable example.
*   **AC:**
    1.  Running `devtools::document()` generates the `.Rd` files with no errors or warnings.
    2.  All exported functions and their parameters are fully documented.
    3.  The `roxygen` examples are self-contained and pass `R CMD check`.
    4.  The `NEWS.md` file is updated to reflect the new feature.