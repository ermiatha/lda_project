

setwd("./data")
filename <- "alzheimer25.sas7bdat"
dat <- read_sas(filename)

# data preparation #############################################################
reshape_dat <- function(dat) {
    dat_long <- dat %>%
        pivot_longer(
            cols = matches("^(bprs|cdrsb|abpet|taupet)\\d+"),
            names_to = c(".value", "time"),
            names_pattern = "(.+)(\\d+)"
        )
    dat_long
}

to_factor <- function(dat) {
    dat <- dat %>%
        mutate(
            sex = as.factor(sex),
            edu = as.factor(edu),
            trial = as.factor(trial),
            job = as.factor(job),
            wzc = as.factor(wzc),
            patid_fct = as.factor(patid)
        )
    dat
}

dat <- to_factor(dat)
dat_long <- reshape_dat(dat)


dat_long <- dat_long %>%
    mutate(
        time = as.factor(time),  # levels 0-6
        time_num = as.numeric(time)   # values 1-7
    )


# baseline per subject (year 0)
base <- dat_long %>%
    filter(time == 0) %>%
    transmute(
        patid, trial,
        # Age = age, Sex = sex, Edu = edu, BMI = bmi, Job = job, Income = inkomen,
        # ADL = adl, WZC = wzc,
        cdrsb0 = cdrsb, 
        abpet0 = abpet, 
        taupet0 = taupet
    )

# add baseline + time to all rows
dat_long <- dat_long %>%
    left_join(base, by = c("patid","trial"))

rm(base)

dat_long <- dat_long %>%
    mutate(
        time_num0 = dat_long$time_num - 1,
        age_bin = cut(
            age,
            breaks = c(-Inf, 60,70, 80, Inf),
            labels = c("<60", "60-70", "70-80", "80+"),
            right = FALSE)
    )


dat_long_scaled <- dat_long
vars_to_scale <- c("time_num", "age", "bmi", "adl", "abpet", "taupet", "bprs")
dat_long_scaled[vars_to_scale] <- sapply(dat_long_scaled[vars_to_scale], scale)

# Define dataset for binary outcome
dat_long_bin <- dat_long %>%
    # Binary outcome: 1 = severe (>10), 0 = non-severe
    mutate(
        cdrsb_bin = as.factor(as.character(ifelse(cdrsb > 10, 1, 0))),
        time_num0 = dat_long$time_num - 1
    )

dat_long_bin_scaled <- dat_long_bin
vars_to_scale <- c("time_num", "age", "bmi", "adl", "abpet", "taupet", "bprs")
dat_long_bin_scaled[vars_to_scale] <- sapply(dat_long_bin_scaled[vars_to_scale], scale)

# Expore Missingness ###########################################################

# Count missing values per variable
missing_summary <- dat %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "Variable", values_to = "N_Missing") %>%
    mutate(
        N_Total = nrow(dat),
        Prop_Missing = N_Missing / N_Total
    ) %>%
    filter(N_Missing > 0) %>%
    arrange(desc(N_Missing))

kable(missing_summary, digits = 3, 
      caption = "Variables with Missing Values") %>%
    kable_styling(bootstrap_options = c("striped", "hover"))

# Missingness by time point for key outcome variables
miss_by_time <- dat_long %>%
    group_by(time_num) %>%
    summarise(
        N_Total = n(),
        BPRS_Missing = sum(is.na(bprs)),
        BPRS_Prop = BPRS_Missing / N_Total,
        CDRSB_Missing = sum(is.na(cdrsb)),
        CDRSB_Prop = CDRSB_Missing / N_Total,
        ABPET_Missing = sum(is.na(abpet)),
        ABPET_Prop = ABPET_Missing / N_Total,
        TAUPET_Missing = sum(is.na(taupet)),
        TAUPET_Prop = TAUPET_Missing / N_Total
    )

kable(miss_by_time, digits = 3,
      caption = "Missingness by Time Point") %>%
    kable_styling(bootstrap_options = c("striped", "hover"))

dat_check <- dat_long %>%
    arrange(patid, time_num) %>%
    group_by(patid) %>%
    mutate(
        miss = is.na(bprs),
        miss_started = cummax(miss))


missingness_pattern <- dat_check %>%
    summarise(
        any_missing = any(miss),
        intermittent = any(!miss & miss_started),
        monotone = any_missing & !intermittent)

missingness_pattern %>% 
    count(monotone)

missingness_pattern %>%
    count(intermittent)

# missingness patterns & missingness per variable
dat %>% select(-patid_fct) %>% 
    md.pattern(rotate.names = T)
dat_long %>% select(-patid_fct, -time_num) %>% 
    md.pattern(rotate.names = T)
sapply(dat_long %>% select(-patid_fct, -time_num), function(x) mean(is.na(x)*100)) 


ggplot(dat_long, aes(x=bprs, fill="density")) + 
    geom_density(alpha=0.5)


# Predict missingness of BPRS/CDRSB with logistic regression

reshape_dat <- function(dat) {
    dat_long <- dat %>%
        pivot_longer(
            cols = matches("^(bprs|cdrsb|abpet|taupet)\\d+"), 
            names_to = c(".value", "time"),
            names_pattern = "(.+)(\\d+)"
        )
    dat_long
}

to_factor <- function(dat_long) {
    dat_long <- dat_long %>%
        mutate(
            time_num = as.numeric(time),
            time = as.factor(time),
            sex = as.factor(sex),
            edu = as.factor(edu),
            trial = as.factor(trial),
            job = as.factor(job),
            wzc = as.factor(wzc)
        )
    dat_long
}

dat_long_base <- reshape_dat(dat)
dat_long_base <- to_factor(dat_long_base)

# separate data frame with only baseline covariates
baseline_covariates <- dat %>%
    select(patid, abpet0, taupet0, bprs0) # Add any other baseline-only variables here

# Join the baseline covariates back to the long data frame
dat_long <- left_join(dat_long_base, baseline_covariates, by = "patid")


# Dropout Model
dat_dropout <- dat_long %>%
    group_by(patid) %>%
    arrange(time_num) %>%
    mutate(
        # Indicator for dropout before the next visit
        dropout_next_visit = ifelse(is.na(lead(bprs)), 1, 0),
        # Lagged outcomes to use as predictors
        BPRS_lag = lag(bprs),
        CDRSB_lag = lag(cdrsb)
    ) %>%
    ungroup() %>%
    # Filter for valid observations to model
    filter(
        !is.na(bprs),         # Must be currently in the study
        time_num < 6,         # Cannot drop out after the last visit
        !is.na(BPRS_lag)      # Must have a previous observation to use as a predictor
    )

# Fit the logistic regression model
dropout_model <- glm(
    dropout_next_visit ~ time_num + age + sex + bmi + wzc + abpet0 + taupet0 + bprs0,
    data = dat_dropout,
    family = binomial(link = "logit")
)

odds_ratios <- tidy(dropout_model, conf.int = TRUE, exponentiate = TRUE)
print(odds_ratios)

# Analysis for Continuous outcome (BPRS) #######################################

# standard LMM
model_lmm <- lmer(bprs ~ time_num0 * (sex + age + adl + job + cdrsb0) + wzc + taupet0 + abpet0 + (1 | trial) + (time_num0 | patid),
                  data = dat_long_scaled, REML = T)  # linear time trend
model_lmm_quadr <- lmer(bprs ~ time_num0 * (sex + age + adl + job + cdrsb0) + I(time_num0^2) + wzc + taupet0 + abpet0 + (1 | trial) + (time_num0 | patid),
                        data = dat_long_scaled, REML = T)  # quadratic time trend
summary(model_lmm, ddf = "lme4")
summary(model_lmm_quadr, ddf = "lme4")

# Set up mice for MI
# prepare dataframe (wide format)
dat_mice <- dat_long_scaled %>%
    select(-patid_fct, -time_num, -time_num0, -age_bin) %>%
    mutate(
        patid = as.integer(patid),  # make sure grouping variables are  integer
        trial = as.integer(trial)
    )
# imputation method
meth <- make.method(dat_mice)
# set method for all incomplete vars except bprs to 0 -> dont impute
meth[c("taupet", "abpet", "cdrsb", "bprs")] <- "2l.pan"
# predictor matrix
pred <- make.predictorMatrix(dat_mice)
pred[, "patid"] <- -2  # 
pred[, "trial"] <- -2
pred[, "taupet0"] <- 0
pred[, "abpet0"]  <- 0
pred[, "cdrsb0"]  <- 0
# make imputations
imp.lmm <- mice(dat_mice, method = meth, predictorMatrix = pred, m = 5, donors = 5, maxit = 10, seed = 1212, print = F) 

# Analyse, pool estimates & save the results
implist <- mids2mitml.list(imp.1)

lmer_fits <- with(
    implist,
    expr = {
        lmerTest::lmer(
            bprs ~ as.numeric(time) * (sex + age + adl + job + cdrsb0) +  # time_num0
                wzc + taupet0 + abpet0 +
                (1 | trial) + (as.numeric(time) | patid),
            REML = T)})

summary(mice::pool(lmer_fits))
results <- testEstimates(lmer_fits, extra.pars = T)
lmer_results_df <- as.data.frame(results$estimates)  # different degrees of freedom obtained than with pool()
save(lmer_results_df, file = "./results/lmer_results.RData") 

# Analysis for Binary outcome (CDRSB) ##########################################
# GLMM:
model_glmm <- glmer(cdrsb_bin ~ time_num0 + I(time_num0^2) + age + sex + bmi + adl + abpet + taupet + 
                        (time_num0 | patid_fct),
                    data = dat_long_bin_scaled,
                    family = binomial(link = "logit"),
                    nAGQ = 1,
                    control = glmerControl(optimizer = "bobyqa"))

summary(model_glmm)

# try out random slopes for quadratic time effects
model_glmm_quadr <- glmer(cdrsb_bin ~ time_num0 + I(time_num0^2) + age + sex + bmi + adl + abpet + taupet + bprs + (time_num0 | patid_fct) + (I(time_num0^2) | patid_fct),
                          data = dat_long_bin_scaled,
                          family = binomial(link = "logit"),
                          nAGQ = 1,
                          control = glmerControl(optimizer = "bobyqa"))

summary(model_glmm_quadr)

# Mice setup
# Convert cdrsb_bin to numeric before imputation
dat_mice <- dat_long_bin_scaled %>%
    select(-patid_fct, -time, -time_num, -age_bin) %>%
    mutate(
        patid = as.integer(patid),
        trial = as.integer(trial),
        cdrsb_bin = as.numeric(as.character(cdrsb_bin))  # Convert factor to 0/1 numeric
    )

# imputation setup with threshold for binary variable
cdrsb_threshold <- 10
meth <- make.method(dat_mice)
meth[c("taupet", "abpet", "cdrsb", "bprs")] <- "2l.pan"
meth["cdrsb_bin"] <- "~as.numeric(cdrsb > 10)"  # Note: > not >=

pred <- make.predictorMatrix(dat_mice)
pred[, ] <- 0
pred[, "patid"] <- -2
pred[, "trial"] <- -2
pred[, c("taupet0", "abpet0", "cdrsb0")] <- 0

# define predictor matrix
pred["cdrsb", c("time_num0", "age", "sex", "edu", "bmi", "bprs")] <- 1
pred["bprs", c("time_num0", "age", "sex", "cdrsb")] <- 1
pred["abpet", c("age", "sex", "edu")] <- 1
pred["taupet", c("age", "sex", "edu")] <- 1
pred["cdrsb_bin", ] <- 0  # no predictors for binary cdrsb


# imputation
imp.glmm <- mice(dat_mice, 
                  method = meth, 
                  predictorMatrix = pred, 
                  m = 5, 
                  maxit = 10, 
                  seed = 1212, 
                  print = FALSE)

# Check that cdrsb_bin was correctly derived
complete(imp.glmm, 1) %>% 
    select(cdrsb, cdrsb_bin) %>% 
    head(20)

library(broom.mixed)  # run this before pool() argument to prevent error message
implist <- mids2mitml.list(imp.glmm)

glmm_fits <- with(
    implist,
    expr = {
        glmer(cdrsb_bin ~ time_num0 + I(time_num0^2) + age + sex + bmi + adl + abpet + taupet + 
                  (time_num0 | patid),  # patid should be factor
              family = binomial(link = "logit"),
              nAGQ = 1,
              control = glmerControl(optimizer = "bobyqa"))})

summary(mice::pool(glmm_fits))
results <- testEstimates(glmm_fits, extra.pars = T)
glmm_results_df <- as.data.frame(results$estimates)  # different degrees of freedom obtained than with pool()


save(glmm_results_df, file = "./results/glmm_results.RData")  # use load() function to load dataframe into environment
length(glmm_fits[[1]]@frame[[1]])  # confirm how many obs were used

# GEE
model_gee <- geeglm(as.numeric(as.character(cdrsb_bin)) ~ time_num0 + I(time_num0^2) + taupet + adl + abpet + age + sex + bmi,
                    id = patid_fct,
                    data = dat_long_bin,
                    family = binomial,
                    corstr = "exchangeable")
summary(model_gee)
exp(coef(model_gee))

implist <- mids2mitml.list(imp.glmm)

gee_fits <- with(
    implist,
    expr = {
        geeglm(as.numeric(as.character(cdrsb_bin)) ~ time_num0 + I(time_num0^2) + taupet + adl + abpet + age + sex + bmi,
               id = patid,
               family = binomial,
               corstr = "exchangeable")
    }
)

summary(mice::pool(gee_fits))
results <- testEstimates(gee_fits, extra.pars = T)
gee_results_df <- as.data.frame(results$estimates) 
round(gee_results_df, 3)

save(gee_results_df, file = "./results/gee_results.RData") 

# Weighted GEE (IPW-GEE)

# Prepare dataset & create weights
dat_ipw <- dat_long_bin %>%
    mutate(
        time = time_num0,
        #time2 = time^2,
        cdrsb_bin_num = as.numeric(as.character(cdrsb_bin)),  # 0/1, NA allowed
        miss_bprs = ifelse(is.na(bprs), 1, 0),
        miss_cdrsb = ifelse(is.na(cdrsb), 1, 0),
        r_miss = ifelse(is.na(cdrsb_bin_num), 0, 1)
    ) %>%
    arrange(patid, time)

dat_ipw <- dat_ipw %>%
    group_by(patid) %>%
    mutate(
        cdrsb_bin_lag = lag(cdrsb_bin_num)
    ) %>%
    ungroup()

miss_glm <- glm(
    r_miss ~ time + age + sex + bmi + adl + abpet0 + taupet0 + cdrsb_bin_lag,
    data = dat_ipw,
    family = binomial
)

dat_ipw$pi <- predict(miss_glm, newdata = dat_ipw, type = "response")

miss_glm_time <- glm(
    r_miss ~ time,
    data = dat_ipw,
    family = binomial
)
dat_ipw$pi_num <- predict(miss_glm_time, newdata = dat_ipw, type = "response")

dat_ipw_compl <- dat_ipw
dat_ipw <- dat_ipw %>%
    mutate(
        pi = pmax(pi, 1e-6),  # ensure prob values are > 0
        pi_num = pmax(pi_num, 1e-6),
        w = pi_num / pi
    ) %>%
    filter(r_miss == 1)   # only observed outcomes

summary(dat_ipw$w)


# estimate gee model with weights
model_wgee <- geeglm(as.numeric(as.character(cdrsb_bin)) ~ time_num0 + I(time_num0^2) + bprs + taupet + adl + abpet + age + sex + bmi,
                     id = patid_fct,
                     data = dat_ipw,
                     family = binomial,
                     corstr = "exchangeable",
                     weights = w)

summary(model_wgee)
nobs(model_wgee)

# expore which covariates predict cdrsb missingness
model_miss_cdrsb <- glm(r_miss ~ time + bprs0 + age + sex + bmi + adl + abpet0 + taupet0 + cdrsb_bin_lag + job + edu + inkomen + wzc,
                        data = dat_ipw_compl,
                        family = binomial(link = "logit"))


# Sensitivity Analysis for MNAR ################################################
### Shift imputed MNAR / delta approach
# predictors & methods
pred <- quickpred(dat_long, exclude = "patid")
meth <- make.method(dat_long)
meth["bprs"] <- "norm"

# arms to shift
shift_arms <- c("5","12","19")

# delta values to try
delta_grid <- c(5, 10, 15, 20)

m <- 10
maxit <- 20

imp0 <- mice(
    dat_long,
    m = 1,
    method = meth,
    predictorMatrix = pred,
    maxit = 0,
    printFlag = FALSE
)

post_template <- imp0$post



results_df <- data.frame(
    delta = numeric(),
    mean_bprs = numeric(),
    sd_bprs = numeric(),
    n_imputed = integer()
)
imp_list <- list()

for (delta in delta_grid) {
    
    message("Running delta = ", delta)
    
    post_k <- post_template
    post_k["bprs"] <- paste0(
        "imp[[j]][ where[[i]], i ] <- imp[[j]][ where[[i]], i ] + ",
        "ifelse( as.character(data$trial[ where[[i]] ]) %in% c(",
        paste0("'", shift_arms, "'", collapse = ","),
        "), ", delta, ", 0 )"
    )
    
    # run MI
    imp_k <- mice(
        dat_long,
        m = m,
        method = meth,
        predictorMatrix = pred,
        post = post_k,
        maxit = maxit,
        printFlag = FALSE
    )
    imp_list[[ as.character(delta) ]] <- imp_k
    
    d1 <- complete(imp_k, 1)
    # imp_k is mids object for one delta scenario
    fits <- with(
        data = imp_k,
        expr = lmer(
            bprs ~ time_num * (sex + age + adl + job + cdrsb) + wzc + taupet + abpet +
                (time_num | trial/patid),
            REML = FALSE
        )
    )
    
    # Convert mira -> list of fitted models
    fitlist <- fits$analyses
    
    # Pool fixed effects (Rubin’s rules)
    pooled <- testEstimates(fitlist, var.comp = FALSE)
    print(pooled)
    
    results_df <- rbind(
        results_df,
        data.frame(
            delta = delta,
            mean_bprs = mean(d1$bprs, na.rm = TRUE),
            sd_bprs   = sd(d1$bprs, na.rm = TRUE),
            n_imputed = sum(imp_k$where[,"bprs"])
        )
    )
}



## Tipping point to pool LMM across imputations
get_mitml_estimates <- function(pooled_obj) {
    # Different mitml versions store the printed table in different slots.
    # We check the common ones.
    candidates <- c("estimates", "results", "table", "coef", "coefficients")
    
    for (nm in candidates) {
        if (!is.null(pooled_obj[[nm]])) {
            tab <- pooled_obj[[nm]]
            tab <- as.data.frame(tab)
            tab$term <- rownames(tab)
            rownames(tab) <- NULL
            
            # Ensure required columns exist (your printed output uses these names)
            if (all(c("Estimate", "Std.Error") %in% names(tab))) {
                return(tab)
            }
        }
    }
    
    # If none found, show structure to help debugging
    stop("Could not extract estimates from pooled_obj. Try str(pooled_obj) to see available slots.")
}

# Pool one delta 
pool_one_delta <- function(mids_obj, delta_value, terms_keep = NULL) {
    
    fits <- suppressWarnings(with(
        data = mids_obj,
        expr = lmer(
            bprs ~ time_num * (sex + age + adl + job + cdrsb) + wzc + taupet + abpet +
                (time_num | trial/patid),
            REML = FALSE,
            control = lmerControl(
                optimizer = "bobyqa",
                optCtrl = list(maxfun = 2e5),
                check.conv.singular = "ignore"
            )
        )
    ))
    
    fitlist <- fits$analyses
    
    pooled <- testEstimates(fitlist, extra.pars = FALSE)
    
    tab <- get_mitml_estimates(pooled) %>%
        mutate(
            delta = delta_value,
            conf.low  = Estimate - 1.96 * Std.Error,
            conf.high = Estimate + 1.96 * Std.Error,
            n_models = length(fitlist),
            n_singular = sum(sapply(fitlist, isSingular, tol = 1e-4))
        )
    
    if (!is.null(terms_keep)) tab <- tab %>% filter(term %in% terms_keep)
    tab
}

# Deltas + term selection 
delta_vals <- sort(as.numeric(names(imp_list)))
first_delta <- as.character(delta_vals[1])

tmp_fit <- suppressWarnings(lmer(
    bprs ~ time_num * (sex + age + adl + job + cdrsb) + wzc + taupet + abpet +
        (time_num | trial/patid),
    data = complete(imp_list[[first_delta]], 1),
    REML = FALSE,
    control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5),
        check.conv.singular = "ignore"
    )
))

coef_names <- names(fixef(tmp_fit))

# Keep time effect + any trial-by-time interactions IF they exist in your parametrization
terms_of_interest <- unique(c(
    "time_num",
    grep("time_num:trial|trial.*:time_num", coef_names, value = TRUE)
))

# If you have no trial*time interactions in the model output, don't facet empty panels:
if (length(terms_of_interest) == 1) {
    message("No trial×time interaction terms found in fixef(). Plotting only time_num.")
}

# Pool across all deltas
pooled_all <- bind_rows(lapply(delta_vals, function(d) {
    pool_one_delta(
        mids_obj = imp_list[[as.character(d)]],
        delta_value = d,
        terms_keep = terms_of_interest
    )
}))


