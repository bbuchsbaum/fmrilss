# Parity with pinned Python GLMsingle (fixtures from
# tools/glmsingle_ref/make_fixtures.py; dev-only, not in the built package).
#
# GLMsingle computes in float32 and forms normal equations, so agreement is
# limited by eps32 * kappa^2 of each voxel's design. Continuous outputs are
# checked against that bound; designs where it exceeds 0.1 are "reference
# unstable" and only reported. Discrete choices must match exactly.

eps32 <- 2^-23

check_glms_scenario <- function(name) {
  fx <- read_glms_fixture(name)
  for (variant in fx$variants) {
    fit <- run_glms_fixture(fx, variant)
    inp <- glms_fixture_inputs(fx)
    g <- function(type, key) glms_ref(fx, variant, type, key)
    kappa <- glms_kappa(fit, inp$extras)[fit$typeb$HRFindex]
    bound <- eps32 * kappa^2
    stable <- bound <= 0.1
    tol_v <- pmax(20 * bound, 1e-5)
    info <- sprintf("%s/%s", name, variant)

    expect_equal(fit$typea$onoffR2, g("typea", "onoffR2"), tolerance = 1e-4, info = info)
    expect_equal(fit$typeb$HRFindex, g("typeb", "HRFindex") + 1, info = info)
    expect_equal(fit$typeb$R2, g("typeb", "R2"), tolerance = 1e-4, info = info)
    for (tp in c("typeb", "typec", "typed")) {
      if (is.null(fit[[tp]])) next
      ours <- unname(fit[[tp]]$betasmd)
      ref <- t(g(tp, "betasmd"))
      expect_identical(is.finite(ours), is.finite(ref), info = paste(info, tp))
      err <- vapply(seq_len(ncol(ours)), function(v) rel_err(ours[, v], ref[, v]), numeric(1))
      err[!is.finite(err)] <- 0
      expect_true(all(err[stable] <= tol_v[stable]),
                  info = sprintf("%s %s: max err/tol %.3g", info, tp, max((err / tol_v)[stable])))
    }
    if (!is.null(fit$typec)) {
      expect_equal(fit$typec$pcnum, as.integer(g("typec", "pcnum")), info = info)
      if (!is.null(fit$typec$xvaltrend)) {
        expect_equal(fit$typec$xvaltrend, g("typec", "xvaltrend"), tolerance = 1e-3, info = info)
      }
    }
    if (!is.null(fit$typed)) {
      ref_frac <- g("typed", "FRACvalue")
      same <- abs(fit$typed$FRACvalue - ref_frac) < 1e-6
      expect_true(all(same[stable]), info = paste(info, "FRACvalue"))
    }
  }
}

for (scenario in c("defaults", "extras", "extras_pc0", "sessions", "unequal",
                   "unrepeated", "fast_events", "collinear_extras", "zero_voxel",
                   "single_frac", "no_library")) {
  local({
    sc <- scenario
    test_that(paste("parity with Python GLMsingle:", sc), {
      skip_if_no_glms_fixtures()
      check_glms_scenario(sc)
    })
  })
}
