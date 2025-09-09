## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)
options(digits = 4)

library(MRTAnalysis)

## ----install_package----------------------------------------------------------
# install.packages("MRTAnalysis")

## ----minimum_example----------------------------------------------------------
set.seed(123)

# Simulate a toy dataset
n <- 20
T_val <- 5
id <- rep(1:n, each = T_val)
dp <- rep(1:T_val, times = n)
A <- rbinom(n * T_val, 1, 0.5)
M <- rbinom(n * T_val, 1, plogis(-0.2 + 0.3 * A + 0.1 * dp))
Y <- ave(0.5 * A + 0.7 * M + 0.2 * dp + rnorm(n * T_val), id) # constant within id

dat <- data.frame(id, dp, A, M, Y)

# Minimal mcee call
fit <- mcee(
    data = dat,
    id = "id", dp = "dp",
    outcome = "Y", treatment = "A", mediator = "M",
    time_varying_effect_form = ~1, # constant-over-time NDEE and NIEE
    control_formula_with_mediator = ~ dp + M, # nuisance adjustment
    control_reg_method = "glm", # default method
    rand_prob = 0.5, # known randomization prob
    verbose = FALSE
)

summary(fit)

## ----data_example-------------------------------------------------------------
data(data_time_varying_mediator_distal_outcome)

dat <- data_time_varying_mediator_distal_outcome

dplyr::glimpse(dat)

dplyr::count(dat, id, name = "Ti") |>
    dplyr::summarise(mean_T = mean(Ti), min_T = min(Ti), max_T = max(Ti))

# Delete some decision points for certain individuals to mimic scenarios
# where not all individuals have the same number of decision points.
dat <- dat[!((dat$id == 1 & dat$dp == 10) | (dat$id == 2 & dat$dp %in% c(9, 10))), ]
dplyr::count(dat, id, name = "Ti") |>
    dplyr::summarise(mean_T = mean(Ti), min_T = min(Ti), max_T = max(Ti))

## ----basid_usage--------------------------------------------------------------
set.seed(1)
fit1 <- mcee(
    data = dat,
    id = "id",
    dp = "dp",
    outcome = "Y",
    treatment = "A",
    mediator = "M",
    availability = "I",
    rand_prob = "p_A",
    time_varying_effect_form = ~1, # NDEE and NIEE are constant over time
    control_formula_with_mediator = ~ dp + M + X, # covariate adjustment
    control_reg_method = "glm", # nuisance learners for q, eta, mu, nu
    verbose = FALSE
)
summary(fit1)

## ----time_varying_effect------------------------------------------------------
fit2 <- mcee(
    data = dat,
    id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
    availability = "I", rand_prob = "p_A",
    time_varying_effect_form = ~dp, # NDEE, NIEE vary linearly in time
    control_formula_with_mediator = ~ dp + M + X,
    control_reg_method = "glm",
    verbose = FALSE
)
summary(fit2) # rows now labeled (Intercept) and dp

## ----other_learners_gam-------------------------------------------------------
# Example: GAM (generalized additive model)
set.seed(2)
fit3 <- mcee(
    data = dat,
    id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
    availability = "I", rand_prob = "p_A",
    time_varying_effect_form = ~dp,
    control_formula_with_mediator = ~ s(dp) + s(M) + s(X), # spline formula for mgcv::gam
    control_reg_method = "gam",
    verbose = FALSE
)
summary(fit3)

## ----specific_dp_only---------------------------------------------------------
fit4 <- mcee(
    data = dat,
    id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
    availability = "I", rand_prob = "p_A",
    time_varying_effect_form = ~1,
    control_formula_with_mediator = ~ dp + M + X,
    control_reg_method = "glm",
    specific_dp_only = c(1, 2),
    verbose = FALSE
)
summary(fit4)

## ----fit1_mcee_fit------------------------------------------------------------
fit1$mcee_fit

## ----lincomb_example1---------------------------------------------------------
# difference between direct and indirect excursion effects
summary(fit1, lincomb_joint = c(1, -1))

## ----lincomb_example2---------------------------------------------------------
fit2 <- mcee(
    data = dat,
    id = "id", dp = "dp", outcome = "Y", treatment = "A", mediator = "M",
    availability = "I", rand_prob = "p_A",
    time_varying_effect_form = ~dp,
    control_formula_with_mediator = ~ dp + M + X,
    control_reg_method = "glm",
    verbose = FALSE
)

## ----lincomb_example2_continued-----------------------------------------------
summary(fit2, lincomb_alpha = c(1, 9), lincomb_beta = c(1, 9))

## ----lincomb_example2_joint---------------------------------------------------
L_joint_t10 <- matrix(c(
    1, 9, # alpha part
    -1, -9 # beta part
), nrow = 1)
summary(fit2, lincomb_joint = L_joint_t10)

## ----inspect_nuisance---------------------------------------------------------
summary(fit1, show_nuisance = TRUE)

head(fit1$nuisance_fitted$mu1)

## -----------------------------------------------------------------------------
# Families (binomial vs Gaussian) are chosen automatically when omitted; for `p` and `q` the default is binomial, for the outcome regressions it is Gaussian.

cfg <- list(
    p   = mcee_config_known("p", dat$p_A), # known randomization prob in MRT
    q   = mcee_config_glm("q", ~ dp + X + M),
    eta = mcee_config_glm("eta", ~ dp + X),
    mu  = mcee_config_glm("mu", ~ dp + X + M),
    nu  = mcee_config_glm("nu", ~ dp + X)
)

fit_gen <- mcee_general(
    data = dat,
    id = "id", dp = "dp", outcome = "Y",
    treatment = "A", mediator = "M",
    availability = "I",
    time_varying_effect_form = ~dp,
    config_p = cfg$p, config_q = cfg$q,
    config_eta = cfg$eta, config_mu = cfg$mu, config_nu = cfg$nu,
    verbose = FALSE
)
summary(fit_gen)

## -----------------------------------------------------------------------------
# Fit nuisance regressions manually: p, q, eta, mu, nu
p1_hat <- dat$p_A # known randomization prob in MRT
p1_hat[dat$I == 0] <- 1 # manually set this to avoid wrapper message
q1_hat <- predict(glm(A ~ dp + X + M, family = binomial(), data = dat[dat$I == 1, ]),
    type = "response", newdata = dat
)
q1_hat[dat$I == 0] <- 1 # manually set this to avoid wrapper message
eta1_hat <- predict(lm(Y ~ dp + X, data = dat[dat$A == dat$I, ]), newdata = dat)
eta0_hat <- predict(lm(Y ~ dp + X, data = dat[dat$A == 0, ]), newdata = dat)
mu1_hat <- predict(lm(Y ~ dp + X + M, data = dat[dat$A == dat$I, ]), newdata = dat)
mu0_hat <- predict(lm(Y ~ dp + X + M, data = dat[dat$A == 0, ]), newdata = dat)
nu1_hat <- predict(lm(mu1 ~ dp + X, data = cbind(dat, mu1 = mu1_hat)[dat$A == 0, ]), newdata = dat)
nu0_hat <- predict(lm(mu0 ~ dp + X, data = cbind(dat, mu0 = mu0_hat)[dat$A == dat$I, ]), newdata = dat)

fit_usr <- mcee_userfit_nuisance(
    data = dat,
    id = "id", dp = "dp", outcome = "Y",
    treatment = "A", mediator = "M",
    availability = "I",
    time_varying_effect_form = ~dp,
    p1 = p1_hat,
    q1 = q1_hat,
    eta1 = eta1_hat, eta0 = eta0_hat,
    mu1 = mu1_hat, mu0 = mu0_hat,
    nu1 = nu1_hat, nu0 = nu0_hat,
    verbose = FALSE
)
summary(fit_usr)

## ----check_equal_1------------------------------------------------------------
fit2$mcee_fit[c("alpha_hat", "beta_hat", "alpha_se", "beta_se")]
fit_gen$mcee_fit[c("alpha_hat", "beta_hat", "alpha_se", "beta_se")]

all.equal(fit2$mcee_fit, fit_gen$mcee_fit, tolerance = 1e-6)

## ----check_equal_2------------------------------------------------------------
all.equal(fit2$mcee_fit, fit_usr$mcee_fit, tolerance = 1e-6)

