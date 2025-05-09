# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

ctmc_n2ll_arma <- function(L, dt, ns, from_to, Xb_q_r, Xb_q_m, p, delta, eq_prec = 1.0e-8, link_r = 1L, a_r = 1.0, l_r = 0.0, u_r = 0.0, link_m = 1L, a_m = 1.0, form = 1L, k = 2.0, norm = TRUE, clip = 0.0) {
    .Call(`_walk_ctmc_n2ll_arma`, L, dt, ns, from_to, Xb_q_r, Xb_q_m, p, delta, eq_prec, link_r, a_r, l_r, u_r, link_m, a_m, form, k, norm, clip)
}

ctmc_predict_arma <- function(L, obs, dt, ns, from_to, Xb_q_r, Xb_q_m, p, delta, eq_prec = 1.0e-8, trunc_tol = 1.0e-8, link_r = 1L, a_r = 1.0, l_r = 0.0, u_r = 0.0, link_m = 1L, a_m = 1.0, form = 1L, k = 2.0, norm = TRUE, clip = 0.0) {
    .Call(`_walk_ctmc_predict_arma`, L, obs, dt, ns, from_to, Xb_q_r, Xb_q_m, p, delta, eq_prec, trunc_tol, link_r, a_r, l_r, u_r, link_m, a_m, form, k, norm, clip)
}

logit <- function(x, L = 0.0, U = 0.0) {
    .Call(`_walk_logit`, x, L, U)
}

soft_plus <- function(x, a = 1.0) {
    .Call(`_walk_soft_plus`, x, a)
}

hard_plus <- function(x) {
    .Call(`_walk_hard_plus`, x)
}

clip_Q <- function(Q, clip) {
    .Call(`_walk_clip_Q`, Q, clip)
}

phi_exp_lnG <- function(phi, lnG, prec = 1.0e-8) {
    .Call(`_walk_phi_exp_lnG`, phi, lnG, prec)
}

load_Q_mult <- function(from_to, Xb_q_r, Xb_q_m, ns, link_r = 1L, a_r = 1.0, l_r = 0.0, u_r = 0.0, link_m = 1L, a_m = 1.0, norm = TRUE, clip = 0.0) {
    .Call(`_walk_load_Q_mult`, from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, l_r, u_r, link_m, a_m, norm, clip)
}

load_Q_add <- function(from_to, Xb_q_r, Xb_q_m, ns, link_r = 1L, a_r = 1.0, link_m = 1L, a_m = 1.0, clip = 0.0) {
    .Call(`_walk_load_Q_add`, from_to, Xb_q_r, Xb_q_m, ns, link_r, a_r, link_m, a_m, clip)
}

load_Q_sde <- function(from_to, Xb_q_r, Xb_q_m, ns, k, a_r = 1.0) {
    .Call(`_walk_load_Q_sde`, from_to, Xb_q_r, Xb_q_m, ns, k, a_r)
}

dense_to_sparse <- function(M, tol = 1.0e-8) {
    .Call(`_walk_dense_to_sparse`, M, tol)
}

my_test <- function(v, Q, prec, renorm = TRUE, t2 = TRUE, checks = TRUE) {
    .Call(`_walk_my_test`, v, Q, prec, renorm, t2, checks)
}

