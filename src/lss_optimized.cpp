#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace arma;

// Robust solve of a small symmetric PSD system (Cholesky, then pseudo-inverse).
static arma::vec small_spd_solve(const arma::mat& M, const arma::vec& b) {
    arma::mat R;
    if (arma::chol(R, M)) {
        arma::vec z = arma::solve(arma::trimatl(R.t()), b);
        return arma::solve(arma::trimatu(R), z);
    }
    return arma::pinv(M) * b;
}

// Build the n x T LSS weight matrix W such that beta = W' Y.
//
// `C` must already be residualized against the confounds. With an empty
// `groups` vector a single pooled "other trials" regressor is used (classic
// LSS); otherwise `groups` holds 1-based trial group codes and one summed
// regressor per group is used (LSS-N), with trial i removed from its group.
// `ridge_x` / `ridge_b` are fractional ridge penalties on the trial and
// other-trial coefficients, scaled by the mean design energy (as in OASIS).
// [[Rcpp::export]]
arma::mat lss_weight_matrix_cpp(const arma::mat& C,
                                const arma::ivec& groups,
                                double eps = 1e-12,
                                double ridge_x = 0.0,
                                double ridge_b = 0.0) {
    const uword n = C.n_rows;
    const uword T = C.n_cols;
    arma::mat W(n, T);

    arma::vec CtC = arma::sum(arma::square(C), 0).t();
    const double lx = ridge_x * arma::mean(CtC);

    if (T == 1) {
        double cc = CtC[0] + lx;
        if (cc <= eps) W.zeros(); else W = C / cc;
        return W;
    }

    if (groups.n_elem == 0) {
        arma::vec total = arma::sum(C, 1);
        double ss_tot = arma::dot(total, total);
        arma::vec CtT = C.t() * total;
        arma::vec bt2 = ss_tot - 2.0 * CtT + CtC;
        arma::vec ctb = CtT - CtC;
        bt2 += ridge_b * arma::mean(bt2);
        bt2.elem(arma::find(bt2 < eps)).fill(arma::datum::inf);
        arma::vec alpha = ctb / bt2;
        arma::vec den = CtC + lx - arma::square(ctb) / bt2;
        den.elem(arma::find(den < eps)).fill(eps);
        arma::vec s = (1.0 + alpha) / den;
        arma::vec u = alpha / den;
        W = C.each_row() % s.t();
        W -= total * u.t();
        return W;
    }

    if (groups.n_elem != T) Rcpp::stop("groups must have one entry per trial");
    const int G = groups.max();
    if (groups.min() < 1) Rcpp::stop("groups must be 1-based codes");
    arma::mat A(n, G, arma::fill::zeros);
    for (uword i = 0; i < T; ++i) A.col(groups[i] - 1) += C.col(i);
    arma::mat AtA = A.t() * A;
    arma::mat AtC = A.t() * C;

    double lb = 0.0;
    if (ridge_b > 0) {
        // Mean diagonal of B_i'B_i over trials and non-empty group columns
        double total = 0.0;
        uword count = 0;
        for (uword i = 0; i < T; ++i) {
            const uword gi = groups[i] - 1;
            for (int g = 0; g < G; ++g) {
                double d = AtA(g, g);
                if ((uword)g == gi) d += CtC[i] - 2.0 * AtC(g, i);
                if (d > eps) { total += d; ++count; }
            }
        }
        if (count > 0) lb = ridge_b * total / count;
    }

    arma::vec s(T, arma::fill::zeros);
    arma::mat U(G, T, arma::fill::zeros);
    for (uword i = 0; i < T; ++i) {
        const uword g = groups[i] - 1;
        arma::vec a_c = AtC.col(i);
        arma::mat Mbb = AtA;
        Mbb.row(g) -= a_c.t();
        Mbb.col(g) -= a_c;
        Mbb(g, g) += CtC[i];
        arma::vec mcb = a_c;
        mcb[g] -= CtC[i];

        arma::uvec keep = arma::find(Mbb.diag() > eps);
        if (keep.n_elem == 0) {
            s[i] = CtC[i] + lx > eps ? 1.0 / (CtC[i] + lx) : 0.0;
            continue;
        }
        arma::vec mk = mcb.elem(keep);
        arma::mat Mk = Mbb.submat(keep, keep);
        Mk.diag() += lb;
        arma::vec h = small_spd_solve(Mk, mk);
        double den = std::max(CtC[i] + lx - arma::dot(mk, h), eps);
        arma::vec hk(G, arma::fill::zeros);
        hk.elem(keep) = h;
        s[i] = (1.0 + hk[g]) / den;
        U.col(i) = -hk / den;
    }
    W = C.each_row() % s.t();
    W += A * U;
    return W;
}

// Blocked W' Y. Voxel blocks are distributed across OpenMP threads, which
// also keeps each Y block cache-resident.
static arma::mat blocked_crossprod(const arma::mat& W, const arma::mat& Y,
                                   uword block_size) {
    const uword V = Y.n_cols;
    arma::mat out(W.n_cols, V);
    const uword n_blocks = (V + block_size - 1) / block_size;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (uword b = 0; b < n_blocks; ++b) {
        const uword j = b * block_size;
        const uword j_end = std::min(j + block_size, V) - 1;
        out.cols(j, j_end) = W.t() * Y.cols(j, j_end);
    }
    return out;
}

//' Fused Single-Pass LSS Solver (C++)
//'
//' Computes Least Squares-Separate (LSS) beta estimates by residualizing the
//' trial design against the confounds, forming the n x T LSS weight matrix
//' once, and applying it to the data with a single matrix product. The data
//' matrix is never residualized: each weight vector already lies in the
//' confound residual space.
//'
//' @param X The confound regressor matrix (n x k).
//' @param Y The data matrix (n x V).
//' @param C The trial-wise design matrix (n x T).
//' @param block_size The number of voxels per OpenMP block when
//'   `use_omp = TRUE`.
//' @param groups Optional 1-based integer trial group codes (LSS-N); NULL for
//'   a single pooled "other trials" regressor.
//' @param use_omp Logical; distribute voxel blocks across OpenMP threads.
//'   Useful with a single-threaded BLAS. With a multithreaded BLAS a single
//'   matrix product is faster.
//' @param ridge_x,ridge_b Fractional ridge penalties on the trial and
//'   other-trial coefficients.
//' @return A T x V matrix of LSS beta estimates.
//' @keywords internal
// [[Rcpp::export]]
arma::mat lss_fused_optim_cpp(const arma::mat& X,
                              const arma::mat& Y,
                              const arma::mat& C,
                              int block_size = 96,
                              SEXP groups = R_NilValue,
                              bool use_omp = true,
                              double ridge_x = 0.0,
                              double ridge_b = 0.0) {
    if (block_size <= 0) {
        Rcpp::stop("block_size must be positive");
    }

    arma::mat C_res = C;
    if (X.n_cols > 0) {
        arma::mat XtX = X.t() * X;
        arma::mat XtXinv;
        bool ok = false;
        try {
            ok = arma::inv_sympd(XtXinv, XtX);
        } catch (...) {
            ok = false;
        }
        if (!ok) XtXinv = arma::pinv(XtX);
        C_res -= X * (XtXinv * (X.t() * C));
    }

    arma::ivec g;
    if (!Rf_isNull(groups)) g = Rcpp::as<arma::ivec>(groups);
    arma::mat W = lss_weight_matrix_cpp(C_res, g, 1e-12, ridge_x, ridge_b);

    if (use_omp) return blocked_crossprod(W, Y, (uword)block_size);
    return W.t() * Y;
}
