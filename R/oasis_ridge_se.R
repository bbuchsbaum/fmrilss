# OASIS ridge + SE helpers
#' @keywords internal
#' @noRd
.oasis_resolve_ridge <- function(pre, ridge_x, ridge_b, ridge_mode = "absolute", K = 1L) {
  ridge_mode <- match.arg(ridge_mode, c("absolute","fractional"))
  
  if (ridge_mode == "absolute") {
    return(list(lx = as.numeric(ridge_x), lb = as.numeric(ridge_b)))
  }
  
  # Fractional mode: scale by mean design energy
  if (K == 1L) {
    # Scale by mean design energy (units: a'a and b'b)
    mx <- mean(pre$d)
    mb <- mean(pre$s)
    return(list(lx = as.numeric(ridge_x) * mx,
                lb = as.numeric(ridge_b) * mb))
  } else {
    # Scale by mean per-basis energy
    DD <- pre$D
    EE <- pre$E
    N  <- dim(DD)[3]
    mx <- mean(vapply(seq_len(N), function(j) mean(diag(DD[,,j])), numeric(1)))
    mb <- mean(vapply(seq_len(N), function(j) mean(diag(EE[,,j])), numeric(1)))
    return(list(lx = as.numeric(ridge_x) * mx,
                lb = as.numeric(ridge_b) * mb))
  }
}

# --- Helper: per-trial SEs from 2x2 normal equations and SSE ---
#' @keywords internal
#' @noRd
.oasis_se_from_norms <- function(d, alpha, s, ridge_x, ridge_b,
                                 RY_norm2, B, G, N_Y, S_Y, dof) {
  # For each trial j and voxel v:
  #   SSE_jv = ||RY||^2 - 2*(β*n1 + γ*n2) + β^2*d + γ^2*e + 2βγ*c
  #   sigma2_jv = SSE_jv / dof
  #   Var(beta) = sigma2_jv * (G^{-1})_{11}, where G = [[d+lambda_x, c],[c, e+lambda_b]]
  N <- nrow(B)
  V <- ncol(B)
  se <- matrix(NA_real_, N, V)

  for (j in seq_len(N)) {
    dj <- d[j] + ridge_x
    ej <- s[j] + ridge_b
    cj <- alpha[j]
    n1 <- N_Y[j, ]
    n2 <- S_Y - n1
    beta  <- B[j, ]
    gamma <- G[j, ]

    SSE <- RY_norm2 - 2*(beta * n1 + gamma * n2) + 
      (beta^2) * d[j] + (gamma^2) * s[j] + 2*beta*gamma*cj
    sigma2 <- pmax(SSE / dof, 0)
    
    # (G^{-1})_{11} = (e) / (d*e - c^2) with ridge included
    denom <- pmax(dj * ej - cj * cj, .Machine$double.eps)
    g11   <- ej / denom
    se[j, ] <- sqrt(sigma2 * g11)
  }
  se
}
