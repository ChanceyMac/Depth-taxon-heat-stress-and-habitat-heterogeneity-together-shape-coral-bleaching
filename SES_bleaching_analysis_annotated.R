## =============================================================================
## SES 2024 CORAL BLEACHING — ANALYSIS CODE (ESSENTIALS, FULLY ANNOTATED)
## Chancey MacDonald & Cilun Djakiman
##
## PURPOSE
##   Stripped, publication-facing version of the analysis. It fits every model
##   reported in the manuscript, produces every main-text figure, generates the
##   supplementary contrast tables, and includes minimum working examples for
##   the ordered-depth (monotonic) sensitivity analysis. Orphan models used
##   during exploration have been removed.
##
## STATISTICAL FRAMEWORK  (matches methods §c of the manuscript)
##   * Family:    Bernoulli with a logit link (any / severe bleaching).
##   * Priors:    brms defaults (flat/uniform on population-level slopes,
##                weakly-informative on the intercept and sd terms).
##   * Chains:    4 chains, 2,000–4,000 iterations, 1,000 warm-up.
##   * Backend:   cmdstanr, with within-chain threading for speed.
##   * Model comparison: LOOIC (leave-one-out cross-validation, elpd + SE).
##   * Random effects: replicate quadrats nested in sites — (1|Site/Rep).
##   * Continuous predictors are z-scored via scale() inside brm() formulas
##     to keep priors on a comparable scale and improve sampler geometry.
##   * Depth is treated as continuous in the main analysis. The monotonic
##     ordered-factor version (mo(Depth_m) on bleach_dat_fac) is the
##     sensitivity analysis.
##
## OUTPUT PLAN
##   Main-text figures:
##       Fig 2  — any/severe bleaching probability at four depths
##       Fig 3  — any/severe bleaching probability across taxa
##       Fig 4  — any/severe bleaching probability, taxa × depth densities
##       Fig 5  — any/severe bleaching along DHW × depth × taxa gradient
##       Fig 6  — any/severe bleaching along HabHet × depth × taxa gradient
##   Supplementary contrast tables:
##       Contr_table_bl_tax_depth.csv         (any, taxa × depth)
##       Contr_table_blsvr_tax_depth.csv      (severe, taxa × depth)
##       Contr_table_dhwXdepth_bl.csv         (any, DHW contrasts within depth)
##       Contr_table_depthXdhw_bl.csv         (any, depth contrasts within DHW)
##       Contr_table_dhwXdepth_blsvr.csv      (severe, DHW contrasts within depth)
##       Contr_table_depthXdhw_blsvr.csv      (severe, depth contrasts within DHW)
##       Contr_table_HabHetXdepth_bl.csv      (any, HabHet contrasts within depth)
##       Contr_table_depthXHabHet_bl.csv      (any, depth contrasts within HabHet)
##       Contr_table_HabHetXdepth_blsvr.csv   (severe, HabHet contrasts within depth)
##       Contr_table_depthXHabHet_blsvr.csv   (severe, depth contrasts within HabHet)
##
## SCRIPT STRUCTURE
##   0.  Setup and data curation
##   1.  Null and coral-cover models         (any + severe)
##   2.  Depth-only models                   (Fig. 2)
##   3.  Taxa-only models                    (Fig. 3)
##   4.  Depth × Taxa models                 (Fig. 4)
##   5.  DHW and Temp overall models
##   6.  Taxa × DHW × Depth 3-way models     (Fig. 5) + coral-cover check
##   7.  HabHet metric selection             (single-covariate models)
##   8.  Taxa × HabHet × Depth 3-way models  (Fig. 6) + coral-cover check
##   9.  Sensitivity: monotonic depth        (minimum examples)
##  10.  Contrast tables for Supp. Materials
## =============================================================================


# ============================================================================
# 0.  SETUP AND DATA CURATION
# ============================================================================

# Load all packages via the companion library file (brms, tidybayes, dplyr,
# ggplot2, ggridges, patchwork, HDInterval, parameters, modelr, marginaleffects,
# paletteer, ...). This is expected to install anything missing.
source("SES_bleaching_library.R")

# ---- Depth colour palette used consistently across every figure ------------ #
# Chosen for accessibility and colour-blindness safety; ordered shallow → deep.
depth_cols <- c("2"  = "#9BC7E0",   # shallow
                "6"  = "#7FB09A",
                "12" = "#6B7A94",
                "18" = "#A88CA8")   # deep

# ---- Load cleaned bleaching data ----------------------------------------------- #
# The raw file contains all in-water surveys. Bleaching_cat is scored on a
# 1–6 scale. 

bleach_dat <- read.csv("bleaching_data_cleaned.csv") %>%
  mutate(Taxa    = factor(Taxa),
         Site    = factor(Site),
         Rep     = factor(Rep))

# ---- Two versions of Depth_m ----------------------------------------------- #
# MAIN ANALYSIS: continuous. Retains order and lets us fit smooth logistic
# curves against a numeric depth axis, which is what the manuscript reports.
bleach_dat <- bleach_dat %>%
  mutate(Depth_m = as.numeric(as.character(Depth_m)))

# SENSITIVITY: ordered factor. Used with brms' mo() monotonic effects
# (Bürkner & Charpentier 2020) to check whether depth relationships were
# meaningfully non-linear.
bleach_dat_fac <- bleach_dat %>%
  mutate(Depth_m = factor(Depth_m,
                          levels = c("2", "6", "12", "18"),
                          ordered = TRUE))

# ---- Overall bleaching prevalence (reported in main text) ------------------ #
# Sanity check of the sample-level rates before modelling.
bleach_dat %>%
  summarise(percent_bleached         = sum(Bleached,     na.rm = TRUE) / n() * 100,
            percent_severly_bleached = sum(Bleached_sev, na.rm = TRUE) / n() * 100)


# ============================================================================
# 1.  NULL MODELS AND CORAL-COVER MODELS
# ============================================================================
# Two jobs here:
#   (a) Establish the overall probability of any / severe bleaching,
#       accounting only for the nested replicate structure — a random
#       intercept for Rep-within-Site absorbs pseudo-replication among
#       colonies scored in the same quadrat.
#   (b) Test whether adding site-level coral cover improves fit (methods §c,
#       para 3). This is what lets the paper say results are 'agnostic to
#       coral cover' — the LOO comparisons below need to be reproducible.

# --- Any bleaching, null model
bl_nullmdl <- brm(Bleached ~ (1|Site/Rep),
                  data = bleach_dat, family = bernoulli,
                  cores = 4, iter = 4000,          # more iters for tail ESS
                  save_pars = save_pars(all = TRUE),
                  backend = "cmdstanr")
# add_criterion attaches LOO. reloo=TRUE refits any observations with high
# Pareto-k so LOO-IC is stable.
bl_nullmdl <- add_criterion(bl_nullmdl, "loo", save_psis = TRUE, reloo = TRUE)

# --- Severe bleaching, null model
blsvr_nullmdl <- brm(Bleached_sev ~ (1|Site/Rep),
                     data = bleach_dat, family = bernoulli,
                     cores = 4,
                     save_pars = save_pars(all = TRUE),
                     backend = "cmdstanr")
blsvr_nullmdl <- add_criterion(blsvr_nullmdl, "loo", save_psis = TRUE, reloo = TRUE)

# --- Any bleaching + coral cover
# Coral cover enters as a single site-level continuous predictor. If LOO shows
# no improvement over bl_nullmdl, we treat cover as non-informative for
# overall bleaching probability.
bl_ccovmdl <- brm(Bleached ~ Coral_cover + (1|Site/Rep),
                  data = bleach_dat, family = bernoulli,
                  cores = 4, iter = 4000,
                  save_pars = save_pars(all = TRUE),
                  backend = "cmdstanr")
bl_ccovmdl <- add_criterion(bl_ccovmdl, "loo", save_psis = TRUE, reloo = TRUE)

# --- Severe bleaching + coral cover
blsvr_covermdl <- brm(Bleached_sev ~ Coral_cover + (1|Site/Rep),
                      data = bleach_dat, family = bernoulli,
                      cores = 4,
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr")
blsvr_covermdl <- add_criterion(blsvr_covermdl, "loo", save_psis = TRUE, reloo = TRUE)

# LOO comparisons: 
loo_compare(bl_nullmdl,    bl_ccovmdl)
loo_compare(blsvr_nullmdl, blsvr_covermdl)


# ============================================================================
# 2.  DEPTH-ONLY MODELS  (Figure 2)
# ============================================================================
# Question: does depth attenuate bleaching probability, when averaged across
# every taxon and site? Continuous depth (2, 6, 12, 18 m); a random-intercept
# structure preserves the nested design.

bl_depthmdl <- brm(Bleached ~ Depth_m + (1|Site/Rep),
                   data = bleach_dat, family = bernoulli,
                   cores = 4,
                   save_pars = save_pars(all = TRUE),
                   backend = "cmdstanr")
bl_depthmdl <- add_criterion(bl_depthmdl, "loo", save_psis = TRUE, reloo = TRUE)

blsvr_depthmdl <- brm(Bleached_sev ~ Depth_m + (1|Site/Rep),
                      data = bleach_dat, family = bernoulli,
                      cores = 4,
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr")
blsvr_depthmdl <- add_criterion(blsvr_depthmdl, "loo", save_psis = TRUE, reloo = TRUE)

loo_compare(bl_nullmdl,    bl_depthmdl)
loo_compare(blsvr_nullmdl, blsvr_depthmdl)
model_parameters(bl_depthmdl,    exponentiate = TRUE, ci = 0.95)
model_parameters(blsvr_depthmdl, exponentiate = TRUE, ci = 0.95)

# ---- Figure 2a: probability distributions of any bleaching by depth -------- #
# add_epred_draws() gives the epred (fitted expected value on the response
# scale) at each of the four survey depths. re_formula=NA excludes the
# random effects so the curve reflects population-level uncertainty only.
fig2a <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18)) %>%
  add_epred_draws(bl_depthmdl, re_formula = NA) %>%
  mutate(Depth_m = factor(Depth_m)) %>%
  ggplot(aes(x = .epred, group = Depth_m, fill = Depth_m)) +
  geom_density(alpha = 0.60) +
  scale_fill_manual(values = depth_cols) +
  theme_tidybayes() +
  labs(x = "Probability of bleaching", fill = "Depth (m)") 

# ---- Figure 2b: probability distributions of severe bleaching by depth ----- #
fig2b <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18)) %>%
  add_epred_draws(blsvr_depthmdl, re_formula = NA) %>%
  mutate(Depth_m = factor(Depth_m)) %>%
  ggplot(aes(x = .epred, group = Depth_m, fill = Depth_m)) +
  geom_density(alpha = 0.60) +
  scale_fill_manual(values = depth_cols) +
  theme_tidybayes() +
  labs(x = "Probability of severe bleaching", fill = "Depth (m)") 

# Combine into a two-panel figure with patchwork
fig2a | fig2b


# ============================================================================
# 3.  TAXA-ONLY MODELS  (Figure 3)
# ============================================================================
# Question: does bleaching probability vary among coral taxa, agnositc to
# depth? adapt_delta and max_treedepth are raised to keep NUTS happy with
# the moderately sparse Seriatopora and Fungia cells.

bl_taxmdl <- brm(Bleached ~ Taxa + (1|Site/Rep),
                 data = bleach_dat, family = bernoulli,
                 cores = 4,
                 control = list(adapt_delta = 0.98, max_treedepth = 16),
                 save_pars = save_pars(all = TRUE),
                 backend = "cmdstanr")
bl_taxmdl <- add_criterion(bl_taxmdl, "loo", save_psis = TRUE, reloo = TRUE)

# Severe: iter = 4000 and thin = 3 keep bulk/tail ESS above 1000 with
# fewer autocorrelated samples in the tails.
blsvr_taxmdl <- brm(Bleached_sev ~ Taxa + (1|Site/Rep),
                    data = bleach_dat, family = bernoulli,
                    cores = 4, iter = 4000, thin = 3,
                    control = list(adapt_delta = 0.98, max_treedepth = 16),
                    save_pars = save_pars(all = TRUE),
                    backend = "cmdstanr")
blsvr_taxmdl <- add_criterion(blsvr_taxmdl, "loo", save_psis = TRUE, reloo = TRUE)

loo_compare(bl_nullmdl,    bl_taxmdl)
loo_compare(blsvr_nullmdl, blsvr_taxmdl)
model_parameters(bl_taxmdl,    exponentiate = TRUE, ci = 0.95)
model_parameters(blsvr_taxmdl, exponentiate = TRUE, ci = 0.95)

# ---- Figure 3a: any bleaching probability by taxa -------------------------- #
fig3a <- bleach_dat %>%
  modelr::data_grid(Taxa = levels(Taxa)) %>%
  add_epred_draws(bl_taxmdl, re_formula = NA) %>%
  ggplot(aes(x = .epred, fill = Taxa)) +
  geom_density(alpha = 0.6) +
  theme_tidybayes() +
  xlim(c(0, 1)) +
  scale_fill_paletteer_d("ggsci::default_jama") +
  labs(x = "Probability of bleaching", y = "Density", fill = "Taxa")

# ---- Figure 3b: severe bleaching probability by taxa ----------------------- #
fig3b <- bleach_dat %>%
  modelr::data_grid(Taxa = levels(Taxa)) %>%
  add_epred_draws(blsvr_taxmdl, re_formula = NA) %>%
  ggplot(aes(x = .epred, fill = Taxa)) +
  geom_density(alpha = 0.6) +
  theme_tidybayes() +
  scale_fill_paletteer_d("ggsci::default_jama") +
  labs(x = "Probability of severe bleaching", y = "Density", fill = "Taxa")

fig3a | fig3b


# ============================================================================
# 4.  TAXA × DEPTH MODELS  (Figure 4)
# ============================================================================
# Question: does the depth attenuation pattern vary among taxa? scale(Depth_m)
# z-scores depth on the fly so the interaction coefficient is on a comparable
# scale to the taxa main effects.

bl_taxxdepthmdl <- brm(Bleached ~ Taxa * scale(Depth_m) + (1|Site/Rep),
                       data = bleach_dat, family = bernoulli,
                       cores = 4,
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr")
# moment_match = FALSE because Pareto-k is fine here without it, and it's much
# faster than the alternative.
bl_taxxdepthmdl <- add_criterion(bl_taxxdepthmdl, "loo",
                                 save_psis = TRUE, reloo = TRUE, moment_match = FALSE)

blsvr_taxxdepthmdl <- brm(Bleached_sev ~ Taxa * scale(Depth_m) + (1|Site/Rep),
                          data = bleach_dat, family = bernoulli,
                          cores = 4,
                          save_pars = save_pars(all = TRUE),
                          backend = "cmdstanr")
blsvr_taxxdepthmdl <- add_criterion(blsvr_taxxdepthmdl, "loo",
                                    save_psis = TRUE, reloo = TRUE, moment_match = FALSE)

# ci = c(0.89, 0.95) returns both intervals: 89% for narrative statements
# in the text, 95% for the figure captions.
model_parameters(bl_taxxdepthmdl,    exponentiate = TRUE, ci = c(0.89, 0.95))
model_parameters(blsvr_taxxdepthmdl, exponentiate = TRUE, ci = c(0.89, 0.95))

# ---- Figure 4a: any bleaching, taxa × depth densities ---------------------- #
# One density strip per taxon-depth cell. Faceting by Taxa (rows) puts each
# taxon on its own row so the reader compares depths within a taxon.
fig4a <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(bl_taxxdepthmdl, re_formula = NA) %>%
  mutate(Depth_m = factor(Depth_m)) %>%
  ggplot(aes(x = .epred, fill = Depth_m)) +
  geom_density(alpha = 0.5) +
  theme_tidybayes() +
  facet_grid(rows = vars(Taxa)) +
  scale_fill_manual(values = depth_cols) +
  labs(x = "Probability of bleaching", y = "Density", fill = "Depth (m)")

# ---- Figure 4b: severe bleaching, taxa × depth densities ------------------- #
fig4b <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(blsvr_taxxdepthmdl, re_formula = NA) %>%
  mutate(Depth_m = factor(Depth_m)) %>%
  ggplot(aes(x = .epred, fill = Depth_m)) +
  geom_density(alpha = 0.5) +
  theme_tidybayes() +
  facet_grid(rows = vars(Taxa)) +
  scale_fill_manual(values = depth_cols) +
  labs(x = "Probability of severe bleaching", y = "Density", fill = "Depth (m)")

fig4a | fig4b


# ============================================================================
# 5.  DHW AND TEMPERATURE OVERALL MODELS
# ============================================================================
# Overall heat-stress effects (methods §c, para 3). Two proxies:
#   * CRW_DHW  — satellite-derived cumulative Degree Heating Weeks
#   * Temp     — in-situ SST at the moment of sampling
# Both are added as single fixed effects to the null-plus-random-effects
# structure so their independent contributions can be quantified.

cor(bleach_dat$CRW_DHW, bleach_dat$Temp, use = "complete.obs")

bl_dhw <- brm(Bleached ~ CRW_DHW + (1|Site/Rep),
              data = bleach_dat, family = bernoulli,
              cores = 4,
              save_pars = save_pars(all = TRUE),
              backend = "cmdstanr")
bl_dhw <- add_criterion(bl_dhw, "loo", save_psis = TRUE, reloo = TRUE)

blsvr_dhw <- brm(Bleached_sev ~ CRW_DHW + (1|Site/Rep),
                 data = bleach_dat, family = bernoulli,
                 cores = 4,
                 save_pars = save_pars(all = TRUE),
                 backend = "cmdstanr")
blsvr_dhw <- add_criterion(blsvr_dhw, "loo", save_psis = TRUE, reloo = TRUE)

# In-situ water temperature at sampling. Correlated with DHW but captures the
# realised local temperature exposure at the moment of scoring.
bl_tempmdl <- brm(Bleached ~ Temp + (1|Site/Rep),
                  data = bleach_dat, family = bernoulli,
                  cores = 4,
                  save_pars = save_pars(all = TRUE),
                  backend = "cmdstanr")
bl_tempmdl <- add_criterion(bl_tempmdl, "loo", save_psis = TRUE, reloo = TRUE)

blsvr_tempmdl <- brm(Bleached_sev ~ Temp + (1|Site/Rep),
                     data = bleach_dat, family = bernoulli,
                     cores = 4,
                     save_pars = save_pars(all = TRUE),
                     backend = "cmdstanr")
blsvr_tempmdl <- add_criterion(blsvr_tempmdl, "loo", save_psis = TRUE, reloo = TRUE)

# How much does fit improve over the null by introducing DHW vs by in-situ temperature?
loo_compare(bl_nullmdl,    bl_dhw,    bl_tempmdl)
loo_compare(blsvr_nullmdl, blsvr_dhw, blsvr_tempmdl)

model_parameters(bl_dhw,        exponentiate = TRUE, ci = 0.95)
model_parameters(blsvr_dhw,     exponentiate = TRUE, ci = 0.95)
model_parameters(bl_tempmdl,    exponentiate = TRUE, ci = 0.95)
model_parameters(blsvr_tempmdl, exponentiate = TRUE, ci = 0.95)


# ============================================================================
# 6.  TAXA × DHW × DEPTH 3-WAY INTERACTION MODELS  (Figure 5)
# ============================================================================
# The core question of the paper: does DHW exposure interact with taxon and
# depth to shape bleaching risk? All continuous predictors are z-scored to
# keep priors on comparable scales; the model needs adapt_delta = 0.98 and
# max_treedepth = 16 for reliable exploration of the interaction ridge.

bl_taxdhwdepthmdl_3x <- brm(Bleached ~ Taxa + scale(CRW_DHW) + scale(Depth_m) +
                              scale(CRW_DHW):Taxa +
                              Taxa:scale(Depth_m) +
                              scale(CRW_DHW):Taxa:scale(Depth_m) +
                              (1|Site/Rep),
                            data = bleach_dat, family = bernoulli,
                            cores = 4,
                            control = list(adapt_delta = 0.98, max_treedepth = 16),
                            save_pars = save_pars(all = TRUE),
                            backend = "cmdstanr")
bl_taxdhwdepthmdl_3x <- add_criterion(bl_taxdhwdepthmdl_3x, "loo",
                                      save_psis = TRUE, reloo = FALSE,
                                      moment_match = FALSE)

blsvr_taxdhwdepthmdl_3x <- brm(Bleached_sev ~ Taxa + scale(CRW_DHW) + scale(Depth_m) +
                                 Taxa:scale(CRW_DHW) +
                                 Taxa:scale(Depth_m) +
                                 Taxa:scale(CRW_DHW):scale(Depth_m) +
                                 (1|Site/Rep),
                               data = bleach_dat, family = bernoulli,
                               cores = 4,
                               control = list(adapt_delta = 0.98, max_treedepth = 16),
                               save_pars = save_pars(all = TRUE),
                               backend = "cmdstanr")
blsvr_taxdhwdepthmdl_3x <- add_criterion(blsvr_taxdhwdepthmdl_3x, "loo",
                                         save_psis = TRUE, reloo = FALSE,
                                         moment_match = FALSE)

# ---- Full diagnostic suite ------------------------------------------------- #
# For methods reporting: posterior-predictive checks (density + ECDF overlay),
# HMC diagnostics (E-BFMI, divergences), and trace/acf/rhat plots.
pp_check(bl_taxdhwdepthmdl_3x,    ndraws = 50)
pp_check(bl_taxdhwdepthmdl_3x,    type = "ecdf_overlay")
pp_check(blsvr_taxdhwdepthmdl_3x, ndraws = 50)
pp_check(blsvr_taxdhwdepthmdl_3x, type = "ecdf_overlay")
rstan::check_hmc_diagnostics(bl_taxdhwdepthmdl_3x$fit)
rstan::check_hmc_diagnostics(blsvr_taxdhwdepthmdl_3x$fit)
mcmc_plot(bl_taxdhwdepthmdl_3x,    type = "trace")
mcmc_plot(bl_taxdhwdepthmdl_3x,    type = "acf")
mcmc_plot(bl_taxdhwdepthmdl_3x,    type = "rhat")
mcmc_plot(blsvr_taxdhwdepthmdl_3x, type = "trace")
mcmc_plot(blsvr_taxdhwdepthmdl_3x, type = "acf")
mcmc_plot(blsvr_taxdhwdepthmdl_3x, type = "rhat")

model_parameters(bl_taxdhwdepthmdl_3x,    exponentiate = TRUE, ci = 0.95)
model_parameters(blsvr_taxdhwdepthmdl_3x, exponentiate = TRUE, ci = 0.95)

# ---- Coral-cover check on the 3-way model (methods §c, para 3) ------------- #
# Adding Coral_cover as an extra fixed effect. If LOO doesn't improve,
# the paper's inference is agnostic to cover — the reported claim.
bl_taxdhwdepthmdl_3x_Cover <- brm(Bleached ~ Taxa + scale(CRW_DHW) + scale(Depth_m) +
                                    Taxa:scale(CRW_DHW) + Taxa:scale(Depth_m) +
                                    Taxa:scale(CRW_DHW):scale(Depth_m) +
                                    Coral_cover + (1|Site/Rep),
                                  data = bleach_dat, family = bernoulli,
                                  cores = 4,
                                  control = list(adapt_delta = 0.98, max_treedepth = 16),
                                  save_pars = save_pars(all = TRUE),
                                  backend = "cmdstanr")
bl_taxdhwdepthmdl_3x_Cover <- add_criterion(bl_taxdhwdepthmdl_3x_Cover, "loo",
                                            save_psis = TRUE, reloo = FALSE,
                                            moment_match = FALSE)

blsvr_taxdhwdepthmdl_3x_Cover <- brm(Bleached_sev ~ Taxa + scale(CRW_DHW) + scale(Depth_m) +
                                       Taxa:scale(CRW_DHW) + Taxa:scale(Depth_m) +
                                       Taxa:scale(CRW_DHW):scale(Depth_m) +
                                       Coral_cover + (1|Site/Rep),
                                     data = bleach_dat, family = bernoulli,
                                     cores = 4,
                                     control = list(adapt_delta = 0.98, max_treedepth = 16),
                                     save_pars = save_pars(all = TRUE),
                                     backend = "cmdstanr")
blsvr_taxdhwdepthmdl_3x_Cover <- add_criterion(blsvr_taxdhwdepthmdl_3x_Cover, "loo",
                                               save_psis = TRUE, reloo = FALSE,
                                               moment_match = FALSE)

loo_compare(bl_taxdhwdepthmdl_3x,    bl_taxdhwdepthmdl_3x_Cover)
loo_compare(blsvr_taxdhwdepthmdl_3x, blsvr_taxdhwdepthmdl_3x_Cover)


# ---- In-situ temperature check on the 3-way model -------------------------- #
# Adding Temp (in-situ at sampling depth) as an extra fixed effect. If LOO
# doesn't improve over the DHW-only 3-way, the reported inference is agnostic
# to instantaneous in-situ temperature — the accumulated heat-stress signal
# (DHW) captures what matters and the moment-of-sampling temperature adds
# nothing beyond it.
bl_taxdhwdepthmdl_3x_Temp <- brm(Bleached ~ Taxa + scale(CRW_DHW) + scale(Depth_m) +
                                   Taxa:scale(CRW_DHW) + Taxa:scale(Depth_m) +
                                   Taxa:scale(CRW_DHW):scale(Depth_m) +
                                   scale(Temp) + (1|Site/Rep),
                                 data = bleach_dat, family = bernoulli,
                                 cores = 4,
                                 control = list(adapt_delta = 0.98, max_treedepth = 16),
                                 save_pars = save_pars(all = TRUE),
                                 backend = "cmdstanr")
blsvr_taxdhwdepthmdl_3x_Temp <- add_criterion(bl_taxdhwdepthmdl_3x_Temp, "loo",
                                           save_psis = TRUE, reloo = FALSE,
                                           moment_match = FALSE)

blsvr_taxdhwdepthmdl_3x_Temp <- brm(Bleached_sev ~ Taxa + scale(CRW_DHW) + scale(Depth_m) +
                                      Taxa:scale(CRW_DHW) + Taxa:scale(Depth_m) +
                                      Taxa:scale(CRW_DHW):scale(Depth_m) +
                                      scale(Temp) + (1|Site/Rep),
                                    data = bleach_dat, family = bernoulli,
                                    cores = 4,
                                    control = list(adapt_delta = 0.98, max_treedepth = 16),
                                    save_pars = save_pars(all = TRUE),
                                    backend = "cmdstanr")
blsvr_taxdhwdepthmdl_3x_Temp <- add_criterion(blsvr_taxdhwdepthmdl_3x_Temp, "loo",
                                              save_psis = TRUE, reloo = FALSE,
                                              moment_match = FALSE)


loo_compare(bl_taxdhwdepthmdl_3x,    )
loo_compare(blsvr_taxdhwdepthmdl_3x, blsvr_taxdhwdepthmdl_3x_Temp)

# ---- Data prep for stacked-point overlays used in Figs 5 & 6 --------------- #
# The Figure 5/6 designs show fitted curves with underlying colony-level
# observations stacked above/below the axis by depth. Each depth gets its
# own horizontal band, with jitter inside the band. This block computes
# the y-coordinate for each colony.
depth_levels <- sort(unique(bleach_dat$Depth_m))
n_depths     <- length(depth_levels)
offset       <- 0.05   # gap between the [0,1] ribbon area and the first band
slot_width   <- 0.05   # thickness of each depth band

# For 'any bleaching': depth_slot 3 (18 m) is furthest from the axis in both
# stacks, so shallow → deep runs top → bottom throughout the whole plot.
bleach_dat_stacked <- bleach_dat %>%
  mutate(depth_slot = match(Depth_m, depth_levels),
         y_stack = if_else(Bleached == 1,
                            1 + offset + (n_depths - depth_slot) * slot_width + slot_width/2,
                           -offset -            (depth_slot - 1) * slot_width - slot_width/2),
         y_stack = y_stack + runif(n(), -slot_width/3, slot_width/3))

# For 'severe bleaching': identical geometry, branching on Bleached_sev.
bleach_dat_stacked_sev <- bleach_dat %>%
  mutate(depth_slot = match(Depth_m, depth_levels),
         y_stack = if_else(Bleached_sev == 1,
                            1 + offset + (n_depths - depth_slot) * slot_width + slot_width/2,
                           -offset -            (depth_slot - 1) * slot_width - slot_width/2),
         y_stack = y_stack + runif(n(), -slot_width/3, slot_width/3))

# ---- Figure 5a: any bleaching along DHW gradient by taxa × depth ----------- #
# Counterfactual curves at four depths, faceted by taxon. Stacked observed
# colony points overlay the fitted ribbon.
fig5a <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    Taxa    = levels(bleach_dat$Taxa),
                    CRW_DHW = seq_range(bleach_dat$CRW_DHW, 100) %>% as.numeric) %>%
  add_epred_draws(bl_taxdhwdepthmdl_3x, re_formula = NA) %>%
  ggplot(aes(y = .epred, x = CRW_DHW, fill = factor(Depth_m))) +
  geom_point(data = bleach_dat_stacked,
             aes(x = CRW_DHW, y = y_stack,
                 colour = factor(Depth_m, levels = depth_levels)),
             alpha = 0.5, size = 0.5,
             position = position_jitter(height = 0, width = 0.005),
             inherit.aes = FALSE, show.legend = FALSE) +
  stat_lineribbon(aes(y = .epred), .width = c(.95), alpha = 0.60, linewidth = 0.3) +
  theme_tidybayes() +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(-0.2, 1.2)) +
  scale_fill_manual(values = depth_cols,   name = "Depth (m)") +
  scale_colour_manual(values = depth_cols, name = "Depth (m)") +
  facet_grid(~Taxa) +
  labs(x = "Degree Heating Weeks", y = "Probability of bleaching", fill = "Depth (m)")

# ---- Figure 5b: severe bleaching along DHW gradient by taxa × depth -------- #
fig5b <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    Taxa    = levels(bleach_dat$Taxa),
                    CRW_DHW = seq_range(bleach_dat$CRW_DHW, 100) %>% as.numeric) %>%
  add_epred_draws(blsvr_taxdhwdepthmdl_3x, re_formula = NA) %>%
  ggplot(aes(y = .epred, x = CRW_DHW, fill = factor(Depth_m))) +
  geom_point(data = bleach_dat_stacked_sev,
             aes(x = CRW_DHW, y = y_stack,
                 colour = factor(Depth_m, levels = depth_levels)),
             alpha = 0.5, size = 0.5,
             position = position_jitter(height = 0, width = 0.005),
             inherit.aes = FALSE, show.legend = FALSE) +
  stat_lineribbon(aes(y = .epred), .width = c(.95), alpha = 0.60, linewidth = 0.3) +
  theme_tidybayes() +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(-0.2, 1.2)) +
  scale_fill_manual(values = depth_cols,   name = "Depth (m)") +
  scale_colour_manual(values = depth_cols, name = "Depth (m)") +
  facet_grid(~Taxa) +
  labs(x = "Degree Heating Weeks", y = "Probability of severe bleaching",
       fill = "Depth (m)")

fig5a / fig5b

# ============================================================================
# 7.  HABITAT HETEROGENEITY (HabHet) METRIC SELECTION
# ============================================================================
# Four candidate metrics were computed from remote-sensed bathymetry:
#   mean_150m, max_150m, mean_1000m, max_1000m
# Each is fit as a single fixed effect against the null-plus-random-effects
# structure. Whichever most improves LOO fit is the metric retained for
# the interaction models. The mean_150m metric was retained (methods §c).

bl_bmean_150 <- brm(Bleached ~ HabHet_mean_150m  + (1|Site/Rep),
                    data = bleach_dat, family = bernoulli,
                    cores = 4, save_pars = save_pars(all = TRUE),
                    backend = "cmdstanr")
bl_bmean_150 <- add_criterion(bl_bmean_150,  "loo", save_psis = TRUE, reloo = TRUE)

bl_bmax_150  <- brm(Bleached ~ HabHet_max_150m   + (1|Site/Rep),
                    data = bleach_dat, family = bernoulli,
                    cores = 4, save_pars = save_pars(all = TRUE),
                    backend = "cmdstanr")
bl_bmax_150  <- add_criterion(bl_bmax_150,   "loo", save_psis = TRUE, reloo = TRUE)

bl_bmean_1000 <- brm(Bleached ~ HabHet_mean_1000m + (1|Site/Rep),
                     data = bleach_dat, family = bernoulli,
                     cores = 4, save_pars = save_pars(all = TRUE),
                     backend = "cmdstanr")
bl_bmean_1000 <- add_criterion(bl_bmean_1000, "loo", save_psis = TRUE, reloo = TRUE)

bl_bmax_1000  <- brm(Bleached ~ HabHet_max_1000m  + (1|Site/Rep),
                     data = bleach_dat, family = bernoulli,
                     cores = 4, save_pars = save_pars(all = TRUE),
                     backend = "cmdstanr")
bl_bmax_1000  <- add_criterion(bl_bmax_1000,  "loo", save_psis = TRUE, reloo = TRUE)

# Pairwise comparisons to the null. The retained metric was bl_bmean_150, the most easily interpreted and spatially similar to our site specific sampling protocol.
loo_compare(bl_nullmdl, bl_bmean_150, bl_bmax_150, bl_bmean_1000, bl_bmax_1000)


# ============================================================================
# 8.  TAXA × HabHet × DEPTH 3-WAY INTERACTION MODELS  (Figure 6)
# ============================================================================
# Same architecture as the DHW 3-way, with HabHet_mean_150m replacing DHW.

bl_taxbmean150depthmdl_3x <- brm(Bleached ~ Taxa + scale(HabHet_mean_150m) + scale(Depth_m) +
                                   Taxa:scale(HabHet_mean_150m) +
                                   Taxa:scale(Depth_m) +
                                   Taxa:scale(HabHet_mean_150m):scale(Depth_m) +
                                   (1|Site/Rep),
                                 data = bleach_dat, family = bernoulli,
                                 cores = 4,
                                 control = list(adapt_delta = 0.98, max_treedepth = 16),
                                 save_pars = save_pars(all = TRUE),
                                 backend = "cmdstanr")
bl_taxbmean150depthmdl_3x <- add_criterion(bl_taxbmean150depthmdl_3x, "loo",
                                           save_psis = TRUE, reloo = FALSE
                                           )

blsvr_taxbmean150depthmdl_3x <- brm(Bleached_sev ~ Taxa + scale(HabHet_mean_150m) + scale(Depth_m) +
                                      Taxa:scale(HabHet_mean_150m) +
                                      Taxa:scale(Depth_m) +
                                      Taxa:scale(HabHet_mean_150m):scale(Depth_m) +
                                      (1|Site/Rep),
                                    data = bleach_dat, family = bernoulli,
                                    cores = 4,
                                    control = list(adapt_delta = 0.98, max_treedepth = 16),
                                    save_pars = save_pars(all = TRUE),
                                    backend = "cmdstanr")
blsvr_taxbmean150depthmdl_3x <- add_criterion(blsvr_taxbmean150depthmdl_3x, "loo",
                                              save_psis = TRUE, reloo = FALSE
                                              )

model_parameters(bl_taxbmean150depthmdl_3x,    exponentiate = TRUE, ci = 0.95)
model_parameters(blsvr_taxbmean150depthmdl_3x, exponentiate = TRUE, ci = 0.95)

# ---- Figure 6a: any bleaching along HabHet gradient by taxa × depth -------- #
fig6a <- bleach_dat %>%
  modelr::data_grid(Depth_m          = c(2, 6, 12, 18),
                    Taxa             = levels(bleach_dat$Taxa),
                    HabHet_mean_150m = seq_range(bleach_dat$HabHet_mean_150m, 100)) %>%
  add_epred_draws(bl_taxbmean150depthmdl_3x, re_formula = NA) %>%
  ggplot(aes(y = .epred, x = HabHet_mean_150m, fill = factor(Depth_m))) +
  geom_point(data = bleach_dat_stacked,
             aes(x = HabHet_mean_150m, y = y_stack,
                 colour = factor(Depth_m, levels = depth_levels)),
             alpha = 0.5, size = 0.5,
             position = position_jitter(height = 0, width = 0.005),
             inherit.aes = FALSE, show.legend = FALSE) +
  stat_lineribbon(aes(y = .epred), .width = c(.95), alpha = 0.60, linewidth = 0.3) +
  theme_tidybayes() +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(-0.2, 1.2)) +
  scale_fill_manual(values = depth_cols,   name = "Depth (m)") +
  scale_colour_manual(values = depth_cols, name = "Depth (m)") +
  facet_grid(~Taxa) +
  labs(x = "Habitat heterogeneity", y = "Probability of bleaching",
       fill = "Depth (m)")

# ---- Figure 6b: severe bleaching along HabHet gradient by taxa × depth ----- #
fig6b <- bleach_dat %>%
  modelr::data_grid(Depth_m          = c(2, 6, 12, 18),
                    Taxa             = levels(bleach_dat$Taxa),
                    HabHet_mean_150m = seq_range(bleach_dat$HabHet_mean_150m, 100)) %>%
  add_epred_draws(blsvr_taxbmean150depthmdl_3x, re_formula = NA) %>%
  ggplot(aes(y = .epred, x = HabHet_mean_150m, fill = factor(Depth_m))) +
  geom_point(data = bleach_dat_stacked_sev,
             aes(x = HabHet_mean_150m, y = y_stack,
                 colour = factor(Depth_m, levels = depth_levels)),
             alpha = 0.5, size = 0.5,
             position = position_jitter(height = 0, width = 0.005),
             inherit.aes = FALSE, show.legend = FALSE) +
  stat_lineribbon(aes(y = .epred), .width = c(.95), alpha = 0.60, linewidth = 0.3) +
  theme_tidybayes() +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(-0.2, 1.2)) +
  scale_fill_manual(values = depth_cols,   name = "Depth (m)") +
  scale_colour_manual(values = depth_cols, name = "Depth (m)") +
  facet_grid(~Taxa) +
  labs(x = "Habitat heterogeneity", y = "Probability of severe bleaching",
       fill = "Depth (m)")

fig6a / fig6b

# ---- Coral-cover check on HabHet 3-way (methods §c, para 3) ---------------- #
# Same logic as the DHW 3-way cover check: if LOO doesn't improve, the
# reported HabHet inferences are agnostic to cover.
bl_taxbmean150depthmdl_3x_Cover <- brm(Bleached ~ Taxa + scale(HabHet_mean_150m) + scale(Depth_m) +
                                         Taxa:scale(HabHet_mean_150m) +
                                         Taxa:scale(Depth_m) +
                                         Taxa:scale(HabHet_mean_150m):scale(Depth_m) +
                                         Coral_cover + (1|Site/Rep),
                                       data = bleach_dat, family = bernoulli,
                                       cores = 4,
                                       control = list(adapt_delta = 0.98, max_treedepth = 16),
                                       save_pars = save_pars(all = TRUE),
                                       backend = "cmdstanr")
bl_taxbmean150depthmdl_3x_Cover <- add_criterion(bl_taxbmean150depthmdl_3x_Cover, "loo",
                                                 save_psis = TRUE, reloo = FALSE
                                                 )

loo_compare(bl_taxbmean150depthmdl_3x, bl_taxbmean150depthmdl_3x_Cover)

blsvr_taxbmean150depthmdl_3x_Cover <- brm(Bleached_sev ~ Taxa + scale(HabHet_mean_150m) + scale(Depth_m) +
                                         Taxa:scale(HabHet_mean_150m) +
                                         Taxa:scale(Depth_m) +
                                         Taxa:scale(HabHet_mean_150m):scale(Depth_m) +
                                         Coral_cover + (1|Site/Rep),
                                       data = bleach_dat, family = bernoulli,
                                       cores = 4,
                                       control = list(adapt_delta = 0.98, max_treedepth = 16),
                                       save_pars = save_pars(all = TRUE),
                                       backend = "cmdstanr")
blsvr_taxbmean150depthmdl_3x_Cover <- add_criterion(blsvr_taxbmean150depthmdl_3x_Cover, "loo",
                                                 save_psis = TRUE, reloo = FALSE
)

loo_compare(blsvr_taxbmean150depthmdl_3x, blsvr_taxbmean150depthmdl_3x_Cover)

# ============================================================================
# 9.  SENSITIVITY ANALYSIS: DEPTH AS ORDERED FACTOR  (minimum examples)
# ============================================================================
# The manuscript reports that treating depth as an ordered factor with mo()
# (Bürkner & Charpentier 2020) improved fit substantially in only one of the
# eight main models (Taxa × DHW × Depth for severe bleaching). Below are
# minimum working examples — depth-only, taxa × depth, and taxa × DHW × depth
# monotonic models — sufficient to reproduce the LOO comparisons in that
# sensitivity paragraph. bleach_dat_fac holds the ordered-factor version of
# Depth_m required by mo().

#NOTE: these models take substanitally longer to run than the previous models. 

# ---- Sensitivity 9.1 — depth-only, monotonic ------------------------------- #
bl_depthmdl_mo <- brm(Bleached ~ mo(Depth_m) + (1|Site/Rep),
                      data = bleach_dat_fac, family = bernoulli,
                      cores = 4, save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr")
bl_depthmdl_mo <- add_criterion(bl_depthmdl_mo, "loo",
                                save_psis = TRUE, reloo = TRUE, moment_match = FALSE)

blsvr_depthmdl_mo <- brm(Bleached_sev ~ mo(Depth_m) + (1|Site/Rep),
                         data = bleach_dat_fac, family = bernoulli,
                         cores = 4, save_pars = save_pars(all = TRUE),
                         backend = "cmdstanr")
blsvr_depthmdl_mo <- add_criterion(blsvr_depthmdl_mo, "loo",
                                   save_psis = TRUE, reloo = TRUE, moment_match = FALSE)

loo_compare(bl_depthmdl,    bl_depthmdl_mo)
loo_compare(blsvr_depthmdl, blsvr_depthmdl_mo)

# ---- Sensitivity 9.2 — taxa × depth, monotonic ----------------------------- #
bl_taxxdepthmdl_mo <- brm(Bleached ~ mo(Depth_m) * Taxa + (1|Site/Rep),
                          data = bleach_dat_fac, family = bernoulli,
                          cores = 4, save_pars = save_pars(all = TRUE),
                          backend = "cmdstanr")
bl_taxxdepthmdl_mo <- add_criterion(bl_taxxdepthmdl_mo, "loo",
                                    save_psis = TRUE, reloo = TRUE, moment_match = FALSE)

blsvr_taxxdepthmdl_mo <- brm(Bleached_sev ~ mo(Depth_m) * Taxa + (1|Site/Rep),
                             data = bleach_dat_fac, family = bernoulli,
                             cores = 4, save_pars = save_pars(all = TRUE),
                             backend = "cmdstanr")
blsvr_taxxdepthmdl_mo <- add_criterion(blsvr_taxxdepthmdl_mo, "loo",
                                       save_psis = TRUE, reloo = TRUE, moment_match = FALSE)

loo_compare(bl_taxxdepthmdl,    bl_taxxdepthmdl_mo)
loo_compare(blsvr_taxxdepthmdl, blsvr_taxxdepthmdl_mo)

# ---- Sensitivity 9.3 — three-way with monotonic depth ---------------------- #
# This is the model that improved substantially for severe bleaching, but
# with the trace/acf issues discussed in the sensitivity paragraph. Diagnostic
# plots run below make that clear when the script is executed.
bl_taxdhwdepthmdl_3x_mo <- brm(Bleached ~ Taxa + scale(CRW_DHW) + mo(Depth_m) +
                                 scale(CRW_DHW):Taxa +
                                 Taxa:mo(Depth_m) +
                                 scale(CRW_DHW):Taxa:mo(Depth_m) +
                                 (1|Site/Rep),
                               data = bleach_dat_fac, family = bernoulli,
                               cores = 4,
                               control = list(adapt_delta = 0.98, max_treedepth = 16),
                               save_pars = save_pars(all = TRUE),
                               backend = "cmdstanr")
bl_taxdhwdepthmdl_3x_mo <- add_criterion(bl_taxdhwdepthmdl_3x_mo, "loo",
                                         save_psis = TRUE, reloo = FALSE, moment_match = FALSE)

blsvr_taxdhwdepthmdl_3x_mo <- brm(Bleached_sev ~ Taxa + scale(CRW_DHW) + mo(Depth_m) +
                                    scale(CRW_DHW):Taxa +
                                    Taxa:mo(Depth_m) +
                                    scale(CRW_DHW):Taxa:mo(Depth_m) +
                                    (1|Site/Rep),
                                  data = bleach_dat_fac, family = bernoulli,
                                  cores = 4,
                                  control = list(adapt_delta = 0.98, max_treedepth = 16),
                                  save_pars = save_pars(all = TRUE),
                                  backend = "cmdstanr")
blsvr_taxdhwdepthmdl_3x_mo <- add_criterion(blsvr_taxdhwdepthmdl_3x_mo, "loo",
                                            save_psis = TRUE, reloo = FALSE, moment_match = FALSE)

loo_compare(bl_taxdhwdepthmdl_3x,    bl_taxdhwdepthmdl_3x_mo)
loo_compare(blsvr_taxdhwdepthmdl_3x, blsvr_taxdhwdepthmdl_3x_mo)

# Trace + acf checks for the severe monotonic 3-way (the noted problem case).
# The methods paragraph on sensitivity references these diagnostics.
mcmc_plot(blsvr_taxdhwdepthmdl_3x_mo, type = "trace")
mcmc_plot(blsvr_taxdhwdepthmdl_3x_mo, type = "acf")
mcmc_plot(blsvr_taxdhwdepthmdl_3x_mo, type = "rhat")


# ============================================================================
# 10.  CONTRAST TABLES FOR SUPPLEMENTARY MATERIALS
# ============================================================================
# Each table is built by pivoting epred draws wide, subtracting to form
# pairwise contrasts, then summarising with median + 89% and 95% HDI.
# These CSVs back the supplementary contrast figures cited in the captions
# of Figs 5 and 6.

# ---- 10.1  Depth × Taxa contrasts (from Fig. 4 models) --------------------- #
contr_table_bl_tax_depth <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(bl_taxxdepthmdl, re_formula = NA) %>%
  ungroup() %>%
  mutate(Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = Depth_m, values_from = .epred, names_prefix = "depth_") %>%
  group_by(Taxa) %>%
  reframe(contrast_2m_6m   = depth_2  - depth_6,
          contrast_2m_12m  = depth_2  - depth_12,
          contrast_2m_18m  = depth_2  - depth_18,
          contrast_6m_12m  = depth_6  - depth_12,
          contrast_6m_18m  = depth_6  - depth_18,
          contrast_12m_18m = depth_12 - depth_18) %>%
  pivot_longer(cols = contrast_2m_6m:contrast_12m_18m,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])
write.csv(contr_table_bl_tax_depth, "Contr_table_bl_tax_depth.csv")

contr_table_blsvr_tax_depth <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(blsvr_taxxdepthmdl, re_formula = NA) %>%
  ungroup() %>%
  mutate(Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = Depth_m, values_from = .epred, names_prefix = "depth_") %>%
  group_by(Taxa) %>%
  reframe(contrast_2m_6m   = depth_2  - depth_6,
          contrast_2m_12m  = depth_2  - depth_12,
          contrast_2m_18m  = depth_2  - depth_18,
          contrast_6m_12m  = depth_6  - depth_12,
          contrast_6m_18m  = depth_6  - depth_18,
          contrast_12m_18m = depth_12 - depth_18) %>%
  pivot_longer(cols = contrast_2m_6m:contrast_12m_18m,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3]) %>% 
  print
write.csv(contr_table_blsvr_tax_depth, "Contr_table_blsvr_tax_depth.csv")


# ---- 10.2  DHW × Depth contrasts from 3-way models (Supp Figs 4–7) --------- #
# 10.2a — DHW contrasts within depth × taxa (any bleaching, Supp Fig 5)
# Levels 1, 4, 6 DHW picked for interpretability: low, moderate, high stress.
contr_table_bl_DHW <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    CRW_DHW = c(1, 4, 6),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(bl_taxdhwdepthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(CRW_DHW = factor(CRW_DHW), Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(CRW_DHW, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = CRW_DHW, values_from = .epred, names_prefix = "dhw_") %>%
  group_by(Taxa, Depth_m) %>%
  reframe(contrast_DHW_1_4 = dhw_1 - dhw_4,
          contrast_DHW_1_6 = dhw_1 - dhw_6,
          contrast_DHW_4_6 = dhw_4 - dhw_6) %>%
  pivot_longer(cols = contrast_DHW_1_4:contrast_DHW_4_6,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, Depth_m, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])%>% 
  print
write.csv(contr_table_bl_DHW, "Contr_table_dhwXdepth_bl.csv")

# 10.2b — Depth contrasts within DHW × taxa (any bleaching, Supp Fig 4)
contr_table_bl_depth <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    CRW_DHW = c(1, 4, 6),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(bl_taxdhwdepthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(CRW_DHW = factor(CRW_DHW), Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(CRW_DHW, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = Depth_m, values_from = .epred, names_prefix = "depth_") %>%
  group_by(Taxa, CRW_DHW) %>%
  reframe(contrast_depth_2_6   = depth_2  - depth_6,
          contrast_depth_2_12  = depth_2  - depth_12,
          contrast_depth_2_18  = depth_2  - depth_18,
          contrast_depth_6_12  = depth_6  - depth_12,
          contrast_depth_6_18  = depth_6  - depth_18,
          contrast_depth_12_18 = depth_12 - depth_18) %>%
  pivot_longer(cols = contrast_depth_2_6:contrast_depth_12_18,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, CRW_DHW, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])
write.csv(contr_table_bl_depth, "Contr_table_depthXdhw_bl.csv")

# 10.2c — DHW contrasts within depth × taxa (severe bleaching, Supp Fig 7)
contr_table_blsvr_DHW <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    CRW_DHW = c(1, 4, 6),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(blsvr_taxdhwdepthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(CRW_DHW = factor(CRW_DHW), Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(CRW_DHW, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = CRW_DHW, values_from = .epred, names_prefix = "dhw_") %>%
  group_by(Taxa, Depth_m) %>%
  reframe(contrast_DHW_1_4 = dhw_1 - dhw_4,
          contrast_DHW_1_6 = dhw_1 - dhw_6,
          contrast_DHW_4_6 = dhw_4 - dhw_6) %>%
  pivot_longer(cols = contrast_DHW_1_4:contrast_DHW_4_6,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, Depth_m, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])
write.csv(contr_table_blsvr_DHW, "Contr_table_dhwXdepth_blsvr.csv")

# 10.2d — Depth contrasts within DHW × taxa (severe bleaching, Supp Fig 6)
contr_table_blsvr_depth <- bleach_dat %>%
  modelr::data_grid(Depth_m = c(2, 6, 12, 18),
                    CRW_DHW = c(1, 4, 6),
                    Taxa    = levels(Taxa)) %>%
  add_epred_draws(blsvr_taxdhwdepthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(CRW_DHW = factor(CRW_DHW), Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(CRW_DHW, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = Depth_m, values_from = .epred, names_prefix = "depth_") %>%
  group_by(Taxa, CRW_DHW) %>%
  reframe(contrast_depth_2_6   = depth_2  - depth_6,
          contrast_depth_2_12  = depth_2  - depth_12,
          contrast_depth_2_18  = depth_2  - depth_18,
          contrast_depth_6_12  = depth_6  - depth_12,
          contrast_depth_6_18  = depth_6  - depth_18,
          contrast_depth_12_18 = depth_12 - depth_18) %>%
  pivot_longer(cols = contrast_depth_2_6:contrast_depth_12_18,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, CRW_DHW, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])
write.csv(contr_table_blsvr_depth, "Contr_table_depthXdhw_blsvr.csv")


# ---- 10.3  HabHet × Depth contrasts from 3-way models (Supp Figs 9–12) ----- #
# Use the observed HabHet quantiles rather than hard-coded values, so the code
# still works if the data update. The 'labels' argument in the mutate below
# assigns human-readable factor levels ("min","q25","q50","q75","max") to the
# five quantile values in order.
hh_q <- quantile(bleach_dat$HabHet_mean_150m)

# 10.3a — HabHet contrasts within depth × taxa (any bleaching, Supp Fig 10)
contr_table_bl_HabHet <- bleach_dat %>%
  modelr::data_grid(Depth_m          = c(2, 6, 12, 18),
                    HabHet_mean_150m = hh_q,
                    Taxa             = levels(Taxa)) %>%
  add_epred_draws(bl_taxbmean150depthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(HabHet_lvl = factor(HabHet_mean_150m,
                             labels = c("min","q25","q50","q75","max")),
         Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(HabHet_lvl, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = HabHet_lvl, values_from = .epred, names_prefix = "hh_") %>%
  group_by(Taxa, Depth_m) %>%
  reframe(contrast_min_q25  = hh_min - hh_q25,
          contrast_min_q50  = hh_min - hh_q50,
          contrast_min_q75  = hh_min - hh_q75,
          contrast_min_max  = hh_min - hh_max,
          contrast_q50_q75  = hh_q50 - hh_q75,
          contrast_q50_max  = hh_q50 - hh_max,
          contrast_q75_max  = hh_q75 - hh_max) %>%
  pivot_longer(cols = contrast_min_q25:contrast_q75_max,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, Depth_m, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])
write.csv(contr_table_bl_HabHet, "Contr_table_HabHetXdepth_bl.csv")

# 10.3b — Depth contrasts within HabHet × taxa (any bleaching, Supp Fig 9)
contr_table_bl_habhet_depth <- bleach_dat %>%
  modelr::data_grid(Depth_m          = c(2, 6, 12, 18),
                    HabHet_mean_150m = hh_q,
                    Taxa             = levels(Taxa)) %>%
  add_epred_draws(bl_taxbmean150depthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(HabHet_lvl = factor(HabHet_mean_150m,
                             labels = c("min","q25","q50","q75","max")),
         Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(HabHet_lvl, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = Depth_m, values_from = .epred, names_prefix = "depth_") %>%
  group_by(Taxa, HabHet_lvl) %>%
  reframe(contrast_depth_2_6   = depth_2  - depth_6,
          contrast_depth_2_12  = depth_2  - depth_12,
          contrast_depth_2_18  = depth_2  - depth_18,
          contrast_depth_6_12  = depth_6  - depth_12,
          contrast_depth_6_18  = depth_6  - depth_18,
          contrast_depth_12_18 = depth_12 - depth_18) %>%
  pivot_longer(cols = contrast_depth_2_6:contrast_depth_12_18,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, HabHet_lvl, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])%>% 
  print
write.csv(contr_table_bl_habhet_depth, "Contr_table_depthXHabHet_bl.csv")

# 10.3c — HabHet contrasts within depth × taxa (severe bleaching, Supp Fig 12)
contr_table_blsvr_HabHet <- bleach_dat %>%
  modelr::data_grid(Depth_m          = c(2, 6, 12, 18),
                    HabHet_mean_150m = hh_q,
                    Taxa             = levels(Taxa)) %>%
  add_epred_draws(blsvr_taxbmean150depthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(HabHet_lvl = factor(HabHet_mean_150m,
                             labels = c("min","q25","q50","q75","max")),
         Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(HabHet_lvl, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = HabHet_lvl, values_from = .epred, names_prefix = "hh_") %>%
  group_by(Taxa, Depth_m) %>%
  reframe(contrast_min_q25  = hh_min - hh_q25,
          contrast_min_q50  = hh_min - hh_q50,
          contrast_min_q75  = hh_min - hh_q75,
          contrast_min_max  = hh_min - hh_max,
          contrast_q50_q75  = hh_q50 - hh_q75,
          contrast_q50_max  = hh_q50 - hh_max,
          contrast_q75_max  = hh_q75 - hh_max) %>%
  pivot_longer(cols = contrast_min_q25:contrast_q75_max,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, Depth_m, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])%>% 
  print
write.csv(contr_table_blsvr_HabHet, "Contr_table_HabHetXdepth_blsvr.csv")

# 10.3d — Depth contrasts within HabHet × taxa (severe bleaching, Supp Fig 11)
contr_table_blsvr_habhet_depth <- bleach_dat %>%
  modelr::data_grid(Depth_m          = c(2, 6, 12, 18),
                    HabHet_mean_150m = hh_q,
                    Taxa             = levels(Taxa)) %>%
  add_epred_draws(blsvr_taxbmean150depthmdl_3x, re_formula = NA) %>%
  ungroup() %>%
  mutate(HabHet_lvl = factor(HabHet_mean_150m,
                             labels = c("min","q25","q50","q75","max")),
         Depth_m = factor(Depth_m), Taxa = factor(Taxa)) %>%
  select(HabHet_lvl, Depth_m, Taxa, .epred, .draw) %>%
  pivot_wider(names_from = Depth_m, values_from = .epred, names_prefix = "depth_") %>%
  group_by(Taxa, HabHet_lvl) %>%
  reframe(contrast_depth_2_6   = depth_2  - depth_6,
          contrast_depth_2_12  = depth_2  - depth_12,
          contrast_depth_2_18  = depth_2  - depth_18,
          contrast_depth_6_12  = depth_6  - depth_12,
          contrast_depth_6_18  = depth_6  - depth_18,
          contrast_depth_12_18 = depth_12 - depth_18) %>%
  pivot_longer(cols = contrast_depth_2_6:contrast_depth_12_18,
               names_to = "Contrast", values_to = ".epred") %>%
  group_by(Taxa, HabHet_lvl, Contrast) %>%
  reframe(Estimate      = quantile(.epred, 0.5),
          HDI_95CI_low  = hdi(.epred, ci = 0.95)[,2],
          HDI_95CI_high = hdi(.epred, ci = 0.95)[,3],
          HDI_89CI_low  = hdi(.epred, ci = 0.89)[,2],
          HDI_89CI_high = hdi(.epred, ci = 0.89)[,3])%>% 
  print
write.csv(contr_table_blsvr_habhet_depth, "Contr_table_depthXHabHet_blsvr.csv")


# ============================================================================
# END OF SCRIPT
# ============================================================================
