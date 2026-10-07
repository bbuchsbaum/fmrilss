// Rank-1 GLM (Pedregosa et al., 2015): one HRF per voxel shared by all
// trials, one amplitude per trial. Fitted by exact alternating least squares
// in the K-dimensional basis space: every quantity needed by both steps is a
// small Gram block of the residualized trial design, so after the single
// product U = X'Y no step touches the n x V data again.
//
// Layout: the residualized trial design X is n x (T K), trial-major with
// basis within trial, so trial i occupies columns i*K .. i*K + K - 1.
// U = X'Y is (T K) x V with the same row blocks.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace arma;

namespace {

// Solve a Gram system, retaining its estimable subspace when singular.
// One common scale preserves the minimum-norm solution and makes rank tests
// invariant to the units of the design. No ridge penalty is introduced.
vec spd_solve(const mat& A, const vec& b, const vec& fallback) {
    const double scale = abs(A.diag()).max();
    if (!(scale > 0.0)) return vec(b.n_elem, fill::zeros);
    const mat As = 0.5 * (A + A.t()) / scale;
    const vec bs = b / scale;
    mat R;
    if (chol(R, As) && min(R.diag()) > 1e-6 * max(R.diag())) {
        vec z = solve(trimatl(R.t()), bs);
        return solve(trimatu(R), z);
    }
    vec values;
    mat vectors;
    if (!eig_sym(values, vectors, As)) return fallback;
    const uvec keep = find(values > 1e-12 * values.max());
    if (keep.is_empty()) return vec(b.n_elem, fill::zeros);
    const mat V = vectors.cols(keep);
    const vec out = V * ((V.t() * bs) / values.elem(keep));
    return out.is_finite() ? out : fallback;
}

}  // namespace

// R1-GLMS: rank-1 GLM with separate (LSS) designs. Trials are partitioned
// into G groups with summed blocks A_g = sum_{j in g} X_j. For trial i in
// group g_i the model is
//   y ~ beta_i X_i h + sum_g r_ig (A_g - [g = g_i] X_i) h,
// i.e. one "other trials" regressor per group (G = 1 is the classic LSS
// design of Mumford et al.); the objective sums the T trial-wise residual
// sums of squares.
//
// Gii: K x K x T, Gii(.,.,i) = X_i'X_i. S: K x K x (T G), slice i*G + g holds
// X_i'A_g. GA: (K G) x (K G) with blocks A_g'A_h. groups: 0-based codes.
// yy = per-voxel residual sum of squares of y. H0: K x V start.
// Returns `other` as a (T G) x V matrix (row i*G + g).
// [[Rcpp::export]]
Rcpp::List r1glms_fit_cpp(const arma::mat& U, const arma::cube& Gii,
                          const arma::cube& S, const arma::mat& GA,
                          const arma::uvec& groups, const arma::vec& yy,
                          const arma::mat& H0, int max_iter = 100,
                          double tol = 1e-7) {
    const uword K = Gii.n_rows;
    const uword T = Gii.n_slices;
    const uword G = GA.n_rows / K;
    const uword V = U.n_cols;
    if (U.n_rows != T * K) Rcpp::stop("U must have T*K rows");
    if (S.n_slices != T * G) Rcpp::stop("S must have T*G slices");
    if (groups.n_elem != T || (T > 0 && groups.max() >= G)) Rcpp::stop("bad groups");
    if (H0.n_rows != K || H0.n_cols != V) Rcpp::stop("H0 must be K x V");
    const double eps = 1e-12;
    const uword KK = K * K;
    const double* gii = Gii.memptr();   // slice i at gii + i*KK
    const double* sp = S.memptr();      // slice i*G+g at sp + (i*G+g)*KK
    // Symmetrized S slices, used by the h step
    std::vector<double> ssym(S.n_elem);
    for (uword sl = 0; sl < S.n_slices; ++sl)
        for (uword c = 0; c < K; ++c)
            for (uword r0 = 0; r0 < K; ++r0)
                ssym[sl * KK + c * K + r0] = sp[sl * KK + c * K + r0] + sp[sl * KK + r0 * K + c];
    // GA blocks stored contiguously: block (p, q) at ga + (p*G + q)*KK
    std::vector<double> ga(G * G * KK);
    for (uword p = 0; p < G; ++p)
        for (uword q = 0; q < G; ++q)
            for (uword c = 0; c < K; ++c)
                for (uword r0 = 0; r0 < K; ++r0)
                    ga[(p * G + q) * KK + c * K + r0] = GA(p * K + r0, q * K + c);

    auto quad = [K](const double* M, const double* x) {
        double out = 0.0;
        for (uword c = 0; c < K; ++c) {
            double col = 0.0;
            for (uword r0 = 0; r0 < K; ++r0) col += M[c * K + r0] * x[r0];
            out += col * x[c];
        }
        return out;
    };

    mat H(K, V), B(T, V), Rm(T * G, V);
    vec obj(V), iters(V);
    ivec conv(V);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
    #endif
    for (uword v = 0; v < V; ++v) {
        const double* u = U.colptr(v);
        vec h = H0.col(v);
        double hn = norm(h);
        if (hn > 0) h /= hn;
        std::vector<double> beta(T, 0.0), r(T * G, 0.0);
        std::vector<double> AY(K * G, 0.0), ay(G), GAh(G * G), sa(G);
        std::vector<double> M((G + 1) * (G + 1)), rhs(G + 1), sol(G + 1);
        std::vector<uword> keep(G + 1);
        std::vector<double> L((G + 1) * (G + 1)), bk(G + 1);
        for (uword i = 0; i < T; ++i)
            for (uword k = 0; k < K; ++k) AY[groups[i] * K + k] += u[i * K + k];
        mat Hm(K, K);
        vec g(K);
        std::vector<double> rr(G * G);

        auto beta_step = [&](const vec& hv) {
            const double* hh = hv.memptr();
            for (uword p = 0; p < G; ++p) {
                double acc = 0.0;
                for (uword k = 0; k < K; ++k) acc += AY[p * K + k] * hh[k];
                ay[p] = acc;
                for (uword q = 0; q < G; ++q) GAh[p * G + q] = quad(&ga[(p * G + q) * KK], hh);
            }
            for (uword i = 0; i < T; ++i) {
                const uword gi = groups[i];
                const double cc = quad(gii + i * KK, hh);
                double cy = 0.0;
                for (uword k = 0; k < K; ++k) cy += u[i * K + k] * hh[k];
                for (uword p = 0; p < G; ++p) sa[p] = quad(sp + (i * G + p) * KK, hh);
                double* ri = &r[i * G];
                if (G == 1) {
                    // Closed-form 2 x 2 solve (classic LSS design)
                    const double cb = sa[0] - cc;
                    const double bb = GAh[0] - 2.0 * sa[0] + cc;
                    const double by = ay[0] - cy;
                    const double det = cc * bb - cb * cb;
                    if (cc > 0.0 && bb > eps * cc && det > eps * cc * bb) {
                        beta[i] = (bb * cy - cb * by) / det;
                        ri[0] = (cc * by - cb * cy) / det;
                    } else {
                        mat gram = {{cc, cb}, {cb, bb}};
                        vec coeff = spd_solve(gram, vec({cy, by}), vec(2, fill::zeros));
                        beta[i] = coeff[0];
                        ri[0] = coeff[1];
                    }
                    continue;
                }
                // Gram of [c, d_1..d_G] with d_g = A_g h - [g = gi] c
                const uword D = G + 1;
                M[0] = cc;
                rhs[0] = cy;
                for (uword p = 0; p < G; ++p) {
                    const double dp = (p == gi) ? 1.0 : 0.0;
                    M[(p + 1) * D] = M[p + 1] = sa[p] - dp * cc;
                    rhs[p + 1] = ay[p] - dp * cy;
                    for (uword q = 0; q < G; ++q) {
                        const double dq = (q == gi) ? 1.0 : 0.0;
                        M[(q + 1) * D + p + 1] =
                            GAh[p * G + q] - dq * sa[p] - dp * sa[q] + dp * dq * cc;
                    }
                }
                // Remove empty columns using a relative, unit-invariant
                // threshold, then try the fast small Cholesky solve.
                double gram_scale = 0.0;
                for (uword p = 0; p < D; ++p) gram_scale = std::max(gram_scale, M[p * D + p]);
                uword nk = 0;
                for (uword p = 0; p < D; ++p)
                    if (M[p * D + p] > eps * gram_scale) keep[nk++] = p;
                // In-place Cholesky (lower, column-major) of the kept block.
                for (uword a2 = 0; a2 < nk; ++a2) {
                    bk[a2] = rhs[keep[a2]];
                    for (uword b2 = 0; b2 < nk; ++b2) L[b2 * nk + a2] = M[keep[b2] * D + keep[a2]];
                }
                bool pd = true;
                for (uword j = 0; j < nk && pd; ++j) {
                    double d = L[j * nk + j];
                    for (uword k2 = 0; k2 < j; ++k2) d -= L[k2 * nk + j] * L[k2 * nk + j];
                    if (!(d > eps * std::abs(L[j * nk + j]))) { pd = false; break; }
                    d = std::sqrt(d);
                    L[j * nk + j] = d;
                    for (uword a2 = j + 1; a2 < nk; ++a2) {
                        double x = L[j * nk + a2];
                        for (uword k2 = 0; k2 < j; ++k2) x -= L[k2 * nk + a2] * L[k2 * nk + j];
                        L[j * nk + a2] = x / d;
                    }
                }
                std::fill(sol.begin(), sol.end(), 0.0);
                if (pd) {
                    for (uword a2 = 0; a2 < nk; ++a2) {          // forward: L z = b
                        double x = bk[a2];
                        for (uword k2 = 0; k2 < a2; ++k2) x -= L[k2 * nk + a2] * bk[k2];
                        bk[a2] = x / L[a2 * nk + a2];
                    }
                    for (uword a2 = nk; a2-- > 0;) {             // backward: L' x = z
                        double x = bk[a2];
                        for (uword k2 = a2 + 1; k2 < nk; ++k2) x -= L[a2 * nk + k2] * bk[k2];
                        bk[a2] = x / L[a2 * nk + a2];
                    }
                    for (uword a2 = 0; a2 < nk; ++a2) sol[keep[a2]] = bk[a2];
                } else {
                    // Dependent nuisance groups must not remove their whole
                    // span: solve for all estimable coefficients together.
                    mat gram(nk, nk);
                    vec target(nk);
                    for (uword a2 = 0; a2 < nk; ++a2) {
                        target[a2] = rhs[keep[a2]];
                        for (uword b2 = 0; b2 < nk; ++b2)
                            gram(a2, b2) = M[keep[b2] * D + keep[a2]];
                    }
                    const vec coeff = spd_solve(gram, target, vec(nk, fill::zeros));
                    for (uword a2 = 0; a2 < nk; ++a2) sol[keep[a2]] = coeff[a2];
                }
                beta[i] = sol[0];
                for (uword p = 0; p < G; ++p) ri[p] = sol[p + 1];
            }
        };

        double prev = datum::inf, cur = datum::inf;
        int it = 0;
        bool converged = false;
        for (it = 1; it <= max_iter; ++it) {
            beta_step(h);
            // h step: A_i = a_i X_i + sum_g r_ig A_g with a_i = beta_i - r_{i,g_i}
            Hm.zeros();
            g.zeros();
            std::fill(rr.begin(), rr.end(), 0.0);
            double* hm = Hm.memptr();
            double* gp = g.memptr();
            for (uword i = 0; i < T; ++i) {
                const double* ri = &r[i * G];
                const double a = beta[i] - ri[groups[i]];
                const double a2 = a * a;
                const double* Gi = gii + i * KK;
                for (uword e = 0; e < KK; ++e) hm[e] += a2 * Gi[e];
                for (uword q = 0; q < G; ++q) {
                    if (ri[q] == 0.0) continue;
                    const double w = a * ri[q];
                    const double* Sq = &ssym[(i * G + q) * KK];
                    for (uword e = 0; e < KK; ++e) hm[e] += w * Sq[e];
                    for (uword p = 0; p < G; ++p) rr[p * G + q] += ri[p] * ri[q];
                    for (uword k = 0; k < K; ++k) gp[k] += ri[q] * AY[q * K + k];
                }
                for (uword k = 0; k < K; ++k) gp[k] += a * u[i * K + k];
            }
            for (uword p = 0; p < G; ++p)
                for (uword q = 0; q < G; ++q) {
                    const double w = rr[p * G + q];
                    if (w == 0.0) continue;
                    const double* Bq = &ga[(p * G + q) * KK];
                    for (uword e = 0; e < KK; ++e) hm[e] += w * Bq[e];
                }
            Hm = 0.5 * (Hm + Hm.t());
            vec h_new = spd_solve(Hm, g, h);
            cur = T * yy[v] - 2.0 * dot(h_new, g) + as_scalar(h_new.t() * Hm * h_new);
            hn = norm(h_new);
            if (!(hn > 0) || !h_new.is_finite()) break;
            h = h_new / hn;  // scale moves into the amplitudes at the next step
            if (std::abs(prev - cur) <= tol * std::max(std::abs(cur), 1e-300)) {
                converged = true;
                break;
            }
            prev = cur;
        }
        beta_step(h);
        // beta_step changes the returned amplitudes after the last h step.
        // Score those final parameters, reusing its h-dependent products.
        cur = T * yy[v];
        for (uword i = 0; i < T; ++i) {
            const double* ri = &r[i * G];
            const double a = beta[i] - ri[groups[i]];
            double cy = 0.0;
            for (uword k = 0; k < K; ++k) cy += u[i * K + k] * h[k];
            cur += a * a * quad(gii + i * KK, h.memptr()) - 2.0 * a * cy;
            for (uword p = 0; p < G; ++p) {
                cur += 2.0 * ri[p] * (a * quad(sp + (i * G + p) * KK, h.memptr()) - ay[p]);
                for (uword q = 0; q < G; ++q) cur += ri[p] * ri[q] * GAh[p * G + q];
            }
        }
        H.col(v) = h;
        for (uword i = 0; i < T; ++i) B(i, v) = beta[i];
        for (uword e = 0; e < T * G; ++e) Rm(e, v) = r[e];
        obj[v] = cur;
        iters[v] = std::min(it, max_iter);
        conv[v] = converged ? 1 : 0;
    }
    return List::create(_["h"] = H, _["beta"] = B, _["other"] = Rm,
                        _["objective"] = obj, _["iterations"] = iters,
                        _["converged"] = conv);
}

// R1-GLM with the joint (LSA) design: y ~ sum_i beta_i X_i h.
// G = X'X is (T K) x (T K). Per voxel and iteration this costs O(T^2 K^2 + T^3),
// so it suits regions of interest rather than whole brains.
// [[Rcpp::export]]
Rcpp::List r1glm_fit_cpp(const arma::mat& U, const arma::mat& G,
                         const arma::vec& yy, const arma::mat& H0,
                         int max_iter = 100, double tol = 1e-7) {
    const uword K = H0.n_rows;
    const uword V = U.n_cols;
    const uword T = U.n_rows / K;
    if (U.n_rows != T * K || G.n_rows != T * K) Rcpp::stop("dimension mismatch");

    mat H(K, V), B(T, V);
    vec obj(V), iters(V);
    ivec conv(V);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
    #endif
    for (uword v = 0; v < V; ++v) {
        vec h = H0.col(v);
        double hn = norm(h);
        if (hn > 0) h /= hn;
        vec beta(T, fill::zeros);
        const vec u = U.col(v);

        auto beta_step = [&](const vec& hh) {
            // D'D(i, j) = h' G_ij h, D'y(i) = h' U_i, with Gh = G (I (x) h)
            mat Gh(T * K, T);
            for (uword j = 0; j < T; ++j) {
                Gh.col(j) = G.cols(j * K, j * K + K - 1) * hh;
            }
            mat DtD(T, T);
            vec Dty(T);
            for (uword i = 0; i < T; ++i) {
                DtD.row(i) = hh.t() * Gh.rows(i * K, i * K + K - 1);
                Dty[i] = dot(hh, u.subvec(i * K, i * K + K - 1));
            }
            DtD = 0.5 * (DtD + DtD.t());
            beta = spd_solve(DtD, Dty, beta);
            return yy[v] - 2.0 * dot(beta, Dty) + as_scalar(beta.t() * DtD * beta);
        };

        double prev = datum::inf, cur = datum::inf;
        int it = 0;
        bool converged = false;
        for (it = 1; it <= max_iter; ++it) {
            beta_step(h);
            // h step: Hm = sum_ij beta_i beta_j G_ij, g = sum_i beta_i U_i
            mat Gb(T * K, K, fill::zeros);
            for (uword j = 0; j < T; ++j) {
                Gb += beta[j] * G.cols(j * K, j * K + K - 1);
            }
            mat Hm(K, K, fill::zeros);
            vec g(K, fill::zeros);
            for (uword i = 0; i < T; ++i) {
                Hm += beta[i] * Gb.rows(i * K, i * K + K - 1);
                g += beta[i] * u.subvec(i * K, i * K + K - 1);
            }
            Hm = 0.5 * (Hm + Hm.t());
            vec h_new = spd_solve(Hm, g, h);
            cur = yy[v] - 2.0 * dot(h_new, g) + as_scalar(h_new.t() * Hm * h_new);
            hn = norm(h_new);
            if (!(hn > 0) || !h_new.is_finite()) break;
            h = h_new / hn;
            if (std::abs(prev - cur) <= tol * std::max(std::abs(cur), 1e-300)) {
                converged = true;
                break;
            }
            prev = cur;
        }
        cur = beta_step(h);
        H.col(v) = h;
        B.col(v) = beta;
        obj[v] = cur;
        iters[v] = std::min(it, max_iter);
        conv[v] = converged ? 1 : 0;
    }
    return List::create(_["h"] = H, _["beta"] = B, _["objective"] = obj,
                        _["iterations"] = iters, _["converged"] = conv);
}

// ---------------------------------------------------------------------------
// Joint quasi-Newton solver for R1-GLMS (the approach of Pedregosa et al.),
// on the same K-space objective as r1glms_fit_cpp, for comparison and as an
// alternative solver. Variables z = [h (K), beta (T), r (T)];
// F(z) = 1/2 sum_i ||y - A_i h||^2 with A_i = (beta_i - r_i) X_i + r_i T_X.
// R's L-BFGS-B (via roptim) is not thread-safe, so voxels run serially.

#include <roptim.h>

namespace {

class R1GLMSObjective : public roptim::Functor {
 public:
    R1GLMSObjective(const arma::cube& Gii, const arma::cube& S, const arma::mat& GTT,
                    const arma::vec& u, double yy)
        : Gii_(Gii), S_(S), GTT_(GTT), u_(u), yy_(yy),
          K_(GTT.n_rows), T_(Gii.n_slices) {
        TY_.zeros(K_);
        for (arma::uword i = 0; i < T_; ++i) TY_ += u_.subvec(i * K_, i * K_ + K_ - 1);
    }

    double operator()(const arma::vec& z) override {
        arma::vec grad;
        return eval(z, grad, false);
    }
    void Gradient(const arma::vec& z, arma::vec& grad) override {
        eval(z, grad, true);
    }

 private:
    double eval(const arma::vec& z, arma::vec& grad, bool want_grad) {
        const arma::vec h = z.head(K_);
        const double tt = arma::as_scalar(h.t() * GTT_ * h);
        const double ty = arma::dot(h, TY_);
        arma::mat Hm(K_, K_, arma::fill::zeros);
        arma::vec g(K_, arma::fill::zeros);
        if (want_grad) grad.zeros(z.n_elem);
        double f = 0.0;
        for (arma::uword i = 0; i < T_; ++i) {
            const double b = z[K_ + i], r = z[K_ + T_ + i], a = b - r;
            const arma::vec ui = u_.subvec(i * K_, i * K_ + K_ - 1);
            const double cc = arma::as_scalar(h.t() * Gii_.slice(i) * h);
            const double ct = arma::as_scalar(h.t() * S_.slice(i) * h);
            const double uy = arma::dot(h, ui);
            f += 0.5 * yy_ - (a * uy + r * ty) +
                 0.5 * (a * a * cc + 2.0 * a * r * ct + r * r * tt);
            if (want_grad) {
                const double da = -uy + a * cc + r * ct;   // dF/da_i
                const double dr = -ty + a * ct + r * tt;   // dF/dr_i at fixed a_i
                grad[K_ + i] = da;
                grad[K_ + T_ + i] = dr - da;
                const arma::mat& Si = S_.slice(i);
                Hm += (a * a) * Gii_.slice(i) + (a * r) * (Si + Si.t()) + (r * r) * GTT_;
                g += a * ui + r * TY_;
            }
        }
        if (want_grad) grad.head(K_) = Hm * h - g;
        return f;
    }

    const arma::cube& Gii_;
    const arma::cube& S_;
    const arma::mat& GTT_;
    const arma::vec u_;
    const double yy_;
    const arma::uword K_, T_;
    arma::vec TY_;
};

}  // namespace

// [[Rcpp::export]]
Rcpp::List r1glms_lbfgs_cpp(const arma::mat& U, const arma::cube& Gii,
                            const arma::cube& S, const arma::mat& GTT,
                            const arma::vec& yy, const arma::mat& H0,
                            const arma::mat& B0, const arma::mat& R0,
                            int maxit = 1000, double factr = 1e7,
                            double pgtol = 0.0) {
    const arma::uword K = GTT.n_rows, T = Gii.n_slices, V = U.n_cols;
    arma::mat H(K, V), B(T, V), Rm(T, V);
    arma::vec obj(V), iters(V);
    arma::ivec conv(V);
    for (arma::uword v = 0; v < V; ++v) {
        R1GLMSObjective fn(Gii, S, GTT, U.col(v), yy[v]);
        arma::vec z = arma::join_cols(H0.col(v), B0.col(v), R0.col(v));
        roptim::Roptim<R1GLMSObjective> opt("L-BFGS-B");
        opt.control.maxit = maxit;
        opt.control.factr = factr;
        opt.control.pgtol = pgtol;
        opt.minimize(fn, z);
        z = opt.par();
        H.col(v) = z.head(K);
        B.col(v) = z.subvec(K, K + T - 1);
        Rm.col(v) = z.tail(T);
        obj[v] = 2.0 * opt.value();
        iters[v] = opt.fncount();
        conv[v] = opt.convergence() == 0 ? 1 : 0;
        if (v % 256 == 0) Rcpp::checkUserInterrupt();
    }
    return Rcpp::List::create(_["h"] = H, _["beta"] = B, _["other"] = Rm,
                              _["objective"] = obj, _["iterations"] = iters,
                              _["converged"] = conv);
}
