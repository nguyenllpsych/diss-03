#####################################
## Helper functions for analyses and eda
## Linh Nguyen                      
## nguyenllpsych@gmail.com         
## November 2023  
#####################################

# Research Question 1: Evidence of Assortative Mating

## H1 Function: Baseline Similarity --------------------------------------------
h1_function <- function(
  # character vector of variables for bivariate correlations 
  var_list, 
  # data frame: each column is a profile, each row is the variable name
  prof_list,
  # data frame name and timepoint (default to baseline correlation)
  # profile correlation is always computed for all waves
  df, time = 1) {
  
  ### BIVARIATE CORRELATIONS ###
  
  # initialize cor table
  cor_tab <- data.frame()
  
  # iterate through variable list
  for(var in var_list) {
    
    # extract vectors from data
    p1 <- df[df$time == time, paste0(var, "_1"), drop = T]
    p2 <- df[df$time == time, paste0(var, "_2"), drop = T]
    
    # extract correlation and p-value
    current_cor <- cor.test(p1, p2, method = "pearson",
                            alternative = "two.sided")
    
    # store sample size for later fisher's z test
    sample_size <- current_cor$parameter + 2
    
    # store values
    cor_tab <- rbind(cor_tab, c(var,
                                sample_size,
                                round(current_cor$estimate, 3), 
                                round(current_cor$p.value, 3),
                                round(current_cor$conf.int[1], 3),
                                round(current_cor$conf.int[2], 3)))
  } # END for var LOOP
  
  # rename cor_tab and ensure variable type
  if(!is.null(var_list)){
    names(cor_tab) <- c("variable", "n", "correlation", "p_value", "LL", "UL")
    cor_tab$correlation <- as.numeric(cor_tab$correlation)
    cor_tab$n <- as.numeric(cor_tab$n)
    cor_tab$p_value <- as.numeric(cor_tab$p_value)
    cor_tab$LL <- as.numeric(cor_tab$LL)
    cor_tab$UL <- as.numeric(cor_tab$UL)
  }
  
  ### PROFILE CORRELATIONS ###
  ### FOR ALL WAVES ###
  
  # extract number of timepoints
  time_list <- sort(unique(df$time))

  # initialize profile cor table
  prof_tab <- data.frame()
  
  # initialize full profile df to be merged back into df
  full_profile_df <- data.frame(
    couple = 1:(nrow(df)/length(unique(df$time)))
  )

  # iterate through profile list
  for(profile in names(prof_list)) {
    
    # variable names for each participant
    p1_list <- paste0(prof_list[, profile], "_1")
    p2_list <- paste0(prof_list[, profile], "_2")
    
    # initialize a profile_df to store profile correlations
    #   ncol = 2 + 3 X # profiles
    #     1 couple identifier
    #     1 time point identifier
    #     3 profile correlations raw/centered/std for each profile
    #   nrow = # couples X # time points
    #   to be merged back into df
    profile_df <- data.frame()
    
    # iterate through time point list
    for (time_point in time_list) {
      
      # extract vectors from data
      p1 <- df %>%
        filter(time == time_point) %>%
        select(all_of(p1_list))
      p2 <- df %>%
        filter(time == time_point) %>%
        select(all_of(p2_list))
      
      # RAW profile correlations
      raw   <- rcorr(t(p1), t(p2), type = "pearson")
      raw_r <- diag(raw$r[c(1:nrow(p1)), c((nrow(p1)+1):(nrow(p1)*2))])
      raw_p <- diag(raw$P[c(1:nrow(p1)), c((nrow(p1)+1):(nrow(p1)*2))])
      
      # GENDER-MEAN-CENTERED profile correlations
      p1_centered <- scale(p1, center = TRUE, scale = FALSE)
      p2_centered <- scale(p2, center = TRUE, scale = FALSE)
      centered   <- rcorr(t(p1_centered), t(p2_centered), type = "pearson")
      centered_r <- diag(centered$r[c(1:nrow(p1)), 
                                    c((nrow(p1)+1):(nrow(p1)*2))])
      centered_p <- diag(centered$P[c(1:nrow(p1)), 
                                    c((nrow(p1)+1):(nrow(p1)*2))])
      
      # STANDARDIZED profile correlations
      p1_std <- scale(p1, center = TRUE, scale = TRUE)
      p2_std <- scale(p2, center = TRUE, scale = TRUE)
      std    <- rcorr(t(p1_std), t(p2_std), type = "pearson")
      std_r <- diag(std$r[c(1:nrow(p1)), c((nrow(p1)+1):(nrow(p1)*2))])
      std_p <- diag(std$P[c(1:nrow(p1)), c((nrow(p1)+1):(nrow(p1)*2))])
      
      # store values
      # proportions of couples with significant profile correlations
      prof_tab <- rbind(prof_tab, c(time_point, profile,
                                    round(mean(raw_p < 0.05, na.rm = T), 5),
                                    round(mean(centered_p < 0.05, na.rm = T), 5),
                                    round(mean(std_p < 0.05, na.rm = T), 5)))
      
      # store values of profile correlations in profile_df
      current_profile_df <- data.frame(
        time = time_point,
        raw_r = raw_r,
        centered_r = centered_r,
        std_r = std_r
      )
      names(current_profile_df)[2] <- paste0(profile, "_raw_r")
      names(current_profile_df)[3] <- paste0(profile, "_centered_r")
      names(current_profile_df)[4] <- paste0(profile, "_std_r")
      
      # append profile_df with previous timepoints
      profile_df <- rbind(profile_df, current_profile_df)
      
    } # END for time_point LOOP
    
    # merge all profile_dfs together
    full_profile_df <- cbind(full_profile_df, profile_df)
  } # END for profile LOOP
  
  # fix full_profile_df names with duplicated time columns
  names(full_profile_df) <- make.names(names(full_profile_df), unique = T)
  full_profile_df <- full_profile_df[, !grepl("time\\.", names(full_profile_df))]

  # rename prof_tab
  if(!is.null(prof_list)) {
    names(prof_tab) <- c("time", "profile", 
                         "raw", "centered", "standardized")
  }
  
  # output
  if(!is.null(prof_list)) {
    return(list(bivariate  = cor_tab,
                profile    = prof_tab,
                profile_df = full_profile_df))
  } else {
    return(list(bivariate  = cor_tab))
  }
  
} # END h1_function

## H2 Function: Difference in Correlations -------------------------------------

h2_function <- function(
  # data.frame with at least 2 columns: bivariate output of h1_function()
  #   variable: variable name
  #   correlation: bivariate correlation between-partner
  cor_tab) {
  
  # extract all variable names
  var_list <- cor_tab$variable
  
  # create data frame with all variable pairs
  cor_compare <- as.data.frame(
    t(combn(var_list, 2))
  )
  cor_compare$V1_cor <- NA
  cor_compare$V2_cor <- NA
  cor_compare$z_stat <- NA
  cor_compare$sig <- NA
  
  # compare each pair
  for(pair in 1:nrow(cor_compare)) {
    
    # grab variable data from cor_tab
    V1 <- cor_tab[cor_tab$variable == cor_compare[pair, "V1"],]
    V2 <- cor_tab[cor_tab$variable == cor_compare[pair, "V2"],]
    
    # fisher's z transformation
    V1z <- fisherz(rho = V1$correlation)
    V2z <- fisherz(rho = V2$correlation)
    
    # solve infinity problem
    if(V1z == Inf){V1z <- 0.99999}
    if(V2z == Inf){V2z <- 0.99999}
    if(V1z == -Inf){V1z <- -0.99999}
    if(V2z == -Inf){V2z <- -0.99999}
      
    # sample sizes
    V1n <- V1$n
    V2n <- V2$n
      
    # z statistic
    z_stat <- (V1z - V2z) / sqrt( (1/(V1n - 3)) + (1/(V2n - 3)) )
    
    # store correlation
    cor_compare[pair, "V1_cor"] <- paste0(V1$correlation, 
                                          " [", V1$LL, " - ", V1$UL,"]")
    cor_compare[pair, "V2_cor"] <- paste0(V2$correlation,
                                          " [", V2$LL, " - ", V2$UL,"]")
    
    # store z test results
    cor_compare[pair, "z_stat"] <- round(z_stat, 3)
    cor_compare[pair, "sig"] <- abs(z_stat) > 1.96
    
  } # END for pair LOOP
  
  # return output
  return(cor_compare)
} # END h2_function

## H3 Function: Longitudinal Similarity ----------------------------------------

#h3_function <- function(var_list, df, baseline = 1) {
#  
#  # create cor_tab dataframe for correlated slopes results
#  cor_tab <- data.frame()
#  
#  # create slopes_tab dataframe for fixed effects estimates
#  slopes_tab <- data.frame()
#  
#  # create slope_df dataframe for extracted slopes from separate gender models
#  slope_df <- data.frame(couple = unique(df$couple),
#                         time = baseline)
#  
#  # loop through variables
#  for (var in var_list) {
#    
#    # var name for each participant
#    p1_var <- paste0(var, "_1")
#    p2_var <- paste0(var, "_2")
#    
#    # fit linear mixed model
#    m1 <- lmer(formula = paste(p1_var, "~ time + (1 + time | couple)"),
#               control = lmerControl(optimizer ="Nelder_Mead"),
#               data = df)
#    m2 <- lmer(formula = paste(p2_var, "~ time + (1 + time | couple)"),
#               control = lmerControl(optimizer ="Nelder_Mead"),
#               data = df)
#    
#    # extract and store fixed effects estimates
#    slopes_tab <- rbind(
#      slopes_tab,
#      c("female", var, round(summary(m1)$coefficients[2,], 3)),
#      c("male", var, round(summary(m2)$coefficients[2,], 3))
#    )
#    
#    # extract fitted slopes
#    slope_1 <- coef(m1)$couple[,2]
#    slope_2 <- coef(m2)$couple[,2]
#    
#    # store fitted slopes
#    slope_df <- merge(slope_df,
#                      data.frame(couple = rownames(coef(m1)$couple),
#                                 slope_1,
#                                 slope_2))
#    names(slope_df)[names(slope_df) == "slope_1"] <- paste0("slope_",
#                                                            var, "_1")
#    names(slope_df)[names(slope_df) == "slope_2"] <- paste0("slope_",
#                                                            var, "_2")
#    
#    # extract correlation and p-value
#    current_cor <- cor.test(slope_1, slope_2, method = "pearson",
#                            alternative = "two.sided")
#    
#    # store values
#    cor_tab <- rbind(cor_tab, c(var, 
#                                current_cor$parameter + 2,
#                                round(current_cor$estimate, 3), 
#                                round(current_cor$p.value, 5),
#                                round(current_cor$conf.int[1], 3),
#                                round(current_cor$conf.int[2], 3)))  
#  } # END for var LOOP
#  
#  # rename cor_tab and ensure variable type
#  names(cor_tab) <- c("variable", "n", "correlation", "p_value", "LL", "UL")
#  cor_tab$correlation <- as.numeric(cor_tab$correlation)
#  cor_tab$p_value <- as.numeric(cor_tab$p_value)
#  cor_tab$LL <- as.numeric(cor_tab$LL)
#  cor_tab$UL <- as.numeric(cor_tab$UL)
#  
#  # rename slopes_tab and ensure variable type
#  names(slopes_tab) <- c("gender", "variable", 
#                         "slope", "SE", "df", "t_value", "p_value")
#  
#  # return cor_tab and slope_df
#  return(list(cor_tab = cor_tab,
#              slopes_tab = slopes_tab,
#              slope_df = slope_df))
#} # END h3_function

h3_function <- function(var_list, df, dir, seed = 202407){
  
  # create an empty data frame to store correlated slopes results
  results_df <- data.frame()
  
  # create and empty data frame to store fixed effects estimates
  slopes_tab <- data.frame()
  
  # create slope_df dataframe for extracted individual random slopes
  slope_df <- data.frame(couple = unique(df$couple))

  # loop through variable list
  for (var in var_list){
    
    # check to see if model output already exists
    file_path <- paste0(dir, "/", var, ".RDS")
    if (file.exists(file_path)){
      
      # grab existing model
      current_mod <- readRDS(file_path)
      
    } else {
      
      # specify models
      M1 <- bf(as.formula(paste0(var, "_1 ~ time + (time | couple)")))
      M2 <- bf(as.formula(paste0(var, "_2 ~ time + (time | couple)")))
      
      # run models
      set.seed(seed)
      current_mod <- brm(M1 + M2 + set_rescor(TRUE), data = df, 
                         cores = 4, chains = 4)
      
      # save model output
      saveRDS(current_mod, file = file_path)

    } # END if else file.exists statement
    
    # grab relevant rows
    current_row <- summary(current_mod)$rescor[, c(
      "Estimate", "Est.Error", "l-95% CI", "u-95% CI")]
    current_slopes <- summary(current_mod)$fixed[3:4, c(
      "Estimate", "Est.Error", "l-95% CI", "u-95% CI")]

    # append results to df
    results_df <- rbind(results_df, current_row)
    slopes_tab <- rbind(
      slopes_tab, cbind(c("female", "male"), current_slopes)
    )
    
    # extract fitted slopes
    col_names <- paste0("Estimate.",
                        gsub("[^[:alnum:]]", "", var), 
                        c("1", "2"), "_time")
      
    slope_1 <- as.data.frame(ranef(current_mod)$couple)[col_names[1]]
    slope_2 <- as.data.frame(ranef(current_mod)$couple)[col_names[2]]
    names(slope_1)[1] <- "slope_1"
    names(slope_2)[1] <- "slope_2"
    
    # store fitted slopes
    slope_df <- merge(slope_df,
                      data.frame(couple = rownames(slope_1),
                                 slope_1,
                                 slope_2))
    names(slope_df)[names(slope_df) == "slope_1"] <- paste0("slope_",
                                                            var, "_1")
    names(slope_df)[names(slope_df) == "slope_2"] <- paste0("slope_",
                                                            var, "_2")

  } # END for var LOOP
  
  # rename results
  names(slopes_tab)[names(slopes_tab) == 'c("female", "male")'] <- "gender"

  return(list(results_df = results_df,
              slopes_tab = slopes_tab,
              slope_df = slope_df))
} # END h3_function

## H4 Function: Perceived/Actual Similarity Comparison -------------------------

h4_function <- function(
  # character vector of all personality variables with self/other reports
  #   naming convention: var_self and var_partner
  perception_list,
  
  # time point (default to baseline) and analytic dataframe
  time = 1, df) {
  
  # create dataframes to store results
  #   similarity_df: raw actual and perceived similarity
  #   compare_df: comparison between actual and perceived similarity
  similarity_df <- compare_df <- data.frame()
  
  # loop through self/other variables
  for (var in perception_list) {
    
    # extract vectors from data
    #   p1self:    female partner's self-perception
    #   p1partner: female partner's perception of their partner
    #   p2self:    male partner's self-perception
    #   p2partner: male partner's perception of their partner
    p1self    <- df[df$time == time, paste0(var, "_self_1"), drop = T]
    p1partner <- df[df$time == time, paste0(var, "_partner_1"), drop = T]
    p2self    <- df[df$time == time, paste0(var, "_self_2"), drop = T]
    p2partner <- df[df$time == time, paste0(var, "_partner_2"), drop = T]
    
    # compute correlation and p-value
    #   actual_sim: between both partners' self-reports
    #   perceived_sim_p1: female partner's self perception and 
    #                     female partner's perception of male partner
    #   perceived_sim_p2: male partner's self perception and 
    #                     male partner's perception of female partner
    actual_sim <- cor.test(p1self, p2self, method = "pearson",
                           alternative = "two.sided")    
    perceived_sim_p1 <- cor.test(p1self, p1partner, method = "pearson",
                                 alternative = "two.sided")    
    perceived_sim_p2 <- cor.test(p2self, p2partner, method = "pearson",
                                 alternative = "two.sided")    
    # store cor, ci, pval
    similarity_df <- rbind(
      similarity_df,
      # actual sim
      c("actual", var, 
        paste0(round(actual_sim$estimate, 3),
               " [", round(actual_sim$conf.int[1], 3),
               " - ", 
               round(actual_sim$conf.int[2], 3) ,"]"),
        round(actual_sim$p.value, 3)),
      c("female-perceived", var,
        paste0(round(perceived_sim_p1$estimate, 3),
               " [", round(perceived_sim_p1$conf.int[1], 3),
               " - ", 
               round(perceived_sim_p1$conf.int[2], 3) ,"]"),
        round(perceived_sim_p1$p.value, 3)),
      c("male-perceived", var,
        paste0(round(perceived_sim_p2$estimate, 3),
               " [", round(perceived_sim_p2$conf.int[1], 3),
               " - ", 
               round(perceived_sim_p2$conf.int[2], 3) ,"]"),
        round(perceived_sim_p2$p.value, 3))
    )
    
    # create data frame with all perception pairs
    sim_list <- as.data.frame(
      matrix(c("actual",           "female-perceived",
               "actual",           "male-perceived",
               "female-perceived", "male-perceived"),
             ncol = 2, byrow = TRUE)
    )
    
    # store current personality variable
    sim_list$personality <- var
    
    # fisher's z transformed of each similarity
    actual_z <- fisherz(rho = actual_sim$estimate)
    female_z <- fisherz(rho = perceived_sim_p1$estimate)
    male_z   <- fisherz(rho = perceived_sim_p2$estimate)
    
    # resolve infinity problem
    if(actual_z == Inf){actual_z <- 0.99999}
    if(female_z == Inf){female_z <- 0.99999}
    if(male_z == Inf){male_z <- 0.99999}
    if(actual_z == -Inf){actual_z <- -0.99999}
    if(female_z == -Inf){female_z <- -0.99999}
    if(male_z == -Inf){male_z <- -0.99999}
    
    # sample sizes
    actual_n <- actual_sim$parameter + 2
    female_n <- perceived_sim_p1$parameter + 2
    male_n   <- perceived_sim_p2$parameter + 2
    
    # z tests
    sim_list[1, "z_stat"] <- (actual_z - female_z) / 
      sqrt( (1/(actual_n - 3)) + (1/(female_n - 3)) )
    sim_list[2, "z_stat"] <- (actual_z - male_z) / 
      sqrt( (1/(actual_n - 3)) + (1/(male_n - 3)) )
    sim_list[3, "z_stat"] <- (female_z - male_z) / 
      sqrt( (1/(female_n - 3)) + (1/(male_n - 3)) )
    
    # significance
    sim_list$sig <- abs(sim_list$z_stat) > 1.96
    
    # rbind current sim_list to compare_df
    compare_df <- rbind(
      compare_df,
      sim_list
    )

  } # END for var LOOP
  
  # return results
  names(similarity_df) <- c("similarity", "personality",
                            "correlation", "p-value")
  return(list(similarity_df = similarity_df,
              compare_df = compare_df))
} # END h4_function

# Research Question 2: Benefits of Assortative Mating --------------------------

## H5+6 Functions: Dyadic Response Surface Analysis ----------------------------

# function for dyadic response surface analysis
# 1. fit full model without gender constrain
# 2. fit gender-constrained model and determine whether it is significantly worse
# 3. extract info for the chosen model
# 4. find support for broad congruence:
#   - negative a4
#   - non significant a3 and a5
#   - allowing main effects: a1 and a2 can be different from 0
drsa_function <- function(
  # character vector of personality variables
  var_list, 
  # character vector of relationship quality variables
  quality_list, 
  # directory to store final model objects,
  dir,
  # data frame for analyses
  df,
  # whether or not to scale the predictors
  # due to small slope values
  scale = FALSE) {
  
  # create dataframe of only baseline data
  df <- df[df$time == 0,]
  
  # initialize empty dataframe to store results
  est_df <- results_df <- data.frame()
  
  # loop through predictor variables
  for(var in var_list) {
    # grab vectors of predictors
    p1_var <- df[, paste0(var, "_1"), drop = T] # female
    p2_var <- df[, paste0(var, "_2"), drop = T] # male
    
    # centering data based on grand mean (Schonbrodt et al., 2018)
    grand_mean <- mean(c(p1_var, p2_var), na.rm = T)
    df$centered_1 <- p1_var - grand_mean
    df$centered_2 <- p2_var - grand_mean
    
    # scale results if too small
    if(scale){
      df$centered_1 <- scale(df$centered_1)
      df$centered_2 <- scale(df$centered_2)
    }
  
    # squaring data
    df$centeredsq_1 <- df$centered_1*df$centered_1
    df$centeredsq_2 <- df$centered_2*df$centered_2
    
    # interaction
    df$centeredint <- df$centered_1*df$centered_2

    # loop through outcome variables
    for(qual in quality_list){

      # check to see if model already exists
      # if yes, grab existing model
      # if no, go through full modeling process
      file_path <- paste0(dir, "/", var, "_", qual, ".RDS")
      if (file.exists(file_path)){
        fit_final_mod <- readRDS(file_path)
      } else {
        
        # full model specification
        full_mod <- paste(
      
          # predicting female quality
          #   b1f: actor; b2f: partner; b3f: actor^2; b4f: int; b5f: partner^2
          paste0(qual, "_1 ~ b1f*centered_1 + b2f*centered_2 +
                 b3f*centeredsq_1 + b4f*centeredint + b5f*centeredsq_2"), "\n",
          
          # predicting male quality
          #   b1m: partner; b2m: actor; b3m: partner^2; b4m: int; b5m: actor^2
          paste0(qual, "_2 ~ b2m*centered_2 + b1m*centered_1 +
                 b5m*centeredsq_2 + b4m*centeredint + b3m*centeredsq_1"), "\n",
          
          # residual correlation
          paste0(qual, "_1 ~~ ", qual, "_2"))
    
        # reduced model specification
        reduced_mod <- paste(
          
          full_mod,
          
          # actor effect equality constraints
          '
          b1f == b2m
          b3f == b5m',
          
          # partner effect equality constraints
          '
          b2f == b1m
          b5f == b3m',
          
          # interaction effect equality constraints
          '
          b4f == b4m'
        )
        
        # fit full mod
        fit_full_mod <- sem(full_mod,
                          data=df, meanstructure=TRUE,
                          estimator="ML", missing="fiml",
                          se="boot", bootstrap =100)
        
        # fit reduced mod
        fit_reduced_mod <- sem(reduced_mod,
                             data=df, meanstructure=TRUE,
                             estimator="ML", missing="fiml",
                             se="boot", bootstrap =100)
        
        # model comparison
        compare_mod <- anova(fit_reduced_mod, fit_full_mod)
        compare_mod_p <- compare_mod$`Pr(>Chisq)`[2]
        
        # add auxiliary params to final chosen model
        if(compare_mod_p < .05) {
          
          # if full model is significantly better
          final_mod <- paste(
            full_mod,
            '
            a1f := b1f + b2f
            a2f := b3f + b4f + b5f
            a3f := b1f - b2f
            a4f := b3f - b4f + b5f
            a5f := b3f - b5f
            
            a1m := b1m + b2m
            a2m := b3m + b4m + b5m
            a3m := b1m - b2m
            a4m := b3m - b4m + b5m
            a5m := b3m - b5m'
          )
        } else {
          
          # if full model is not significantly better
          final_mod <- paste(
            reduced_mod,
            '
            a1f := b1f + b2f
            a2f := b3f + b4f + b5f
            a3f := b1f - b2f
            a4f := b3f - b4f + b5f
            a5f := b3f - b5f
            '
          )
        }
        
        # fit final model
        fit_final_mod <- sem(final_mod,
                             data=df, meanstructure=TRUE,
                             estimator="ML", missing="fiml",
                             se="boot", bootstrap=100)
        
        # save final model
        saveRDS(fit_final_mod, file = paste0(file_path))
      } # END else file_path exists STATEMENT

      # extract summary
      est <- as.data.frame(unclass(standardizedSolution(fit_final_mod))) %>% 
        filter(op == ":=") %>%
        mutate(predictor = var,
               outcome = qual) %>%
        select(predictor, outcome, label, 
               est.std, se, pvalue, ci.lower, ci.upper)
      
      # append estimates to est_df
      est_df <- rbind(est_df,
                      est)
      
      # extract congruence results
      actor_partner_f <-
        est[est$label=="a1f", "pvalue"] < .05 | 
        est[est$label=="a2f", "pvalue"] < .05
      broad_congruence_f <- 
        # negative and significant a4
        est[est$label == "a4f", "est.std"] < 0 & 
          est[est$label == "a4f", "pvalue"] < .05 &
        # nonsignificant a3 and a5
        est[est$label == "a3f", "pvalue"] > .05 &
          est[est$label == "a5f", "pvalue"] > .05
      strict_congruence_f <-
        broad_congruence_f & !actor_partner_f
      
      # store congruence results in results_df
      if(nrow(est) == 10){
        # if full model is used -> store both gender results
        actor_partner_m <-
          est[est$label=="a1m", "pvalue"] < .05 | 
          est[est$label=="a2m", "pvalue"] < .05
        broad_congruence_m <- 
          # negative and significant a4
          est[est$label == "a4m", "est.std"] < 0 & 
          est[est$label == "a4m", "pvalue"] < .05 &
          # nonsignificant a3 and a5
          est[est$label == "a3m", "pvalue"] > .05 &
          est[est$label == "a5m", "pvalue"] > .05
        strict_congruence_m <-
          broad_congruence_m & !actor_partner_m

        results <- data.frame(
          predictor = var,
          outcome = qual,
          model = "full",
          sex = c("female", "male"),
          actor_partner = c(actor_partner_f, actor_partner_m),
          broad_congruence = c(broad_congruence_f, broad_congruence_m),
          strict_congruence = c(strict_congruence_f, strict_congruence_m)
        )
      } else {
        # if reduced model is used -> store only f params (they're the same)
        results <- data.frame(
          predictor = var,
          outcome = qual,
          model = "reduced",
          sex = NA,
          actor_partner = actor_partner_f,
          broad_congruence = broad_congruence_f,
          strict_congruence = strict_congruence_f
        )
      } # END if else STATEMENT
      results_df <- rbind(results_df, results)
      
    } # END for qual LOOP
  } # END for var LOOP

  # return output
  return(list(est_df = est_df,
              results_df = results_df))
} #END drsa_function 

# Miscellaneous Functions ------------------------------------------------------

# function to test significance, 2-tailed paired t-test
# return 0/1 for signif
signif <- function(var, var_sex, data, time = 0, female_male) {
  pval <- t.test(x = data[which(data$time == time & data[var_sex] == female_male[1]), 
                                 var],
                 y = data[which(data$time == time & data[var_sex] == female_male[2]), 
                                 var],
                 paired = TRUE,
                 alternative = "two.sided",
                 conf.level = 0.95)$p.value
  signif <- isTRUE(pval < 0.05)
  return(signif)
} # END signif

# function to test significant change, linear mixed model slope
# return 0/1 for signif
signifc <- function(var, data) {
  formula <- as.formula(paste(var, "~ time", "+ (1 + time | IDg)"))
  mod <- summary(lmer(formula = formula,
                      control = lmerControl(optimizer ="Nelder_Mead"),
                      data = data))
  
  signif <- isTRUE(mod$coefficients[2,5] < 0.05)
  return(signif)
} # END signifc

# function for histogram
plot_hist <- function(var, var_name, data, time = 0, bin_width,
                      var_sex, female_male) {
  # custom colors for female and males
  colors <- c("#E69F00", "#009E73")
  names(colors) <- c("female", "male")
  color_text <- glue::glue(
    'for <span style = color:{colors["female"]}>**female scores**</span> ',
    'and ',
    '<span style = color:{colors["male"]}>**male scores**</span>'
  )
  
  p <- ggplot(data = data[which(data$time == time), ], 
       aes(!!sym(var))) +
    geom_histogram(binwidth = bin_width, fill = "#710c0c") + 
    geom_vline(xintercept = 
                 mean(data[which(data$time == 0 & data[var_sex] == female_male[1]), var],
                      na.rm = TRUE),
               color = colors["female"], size = 1) +
    geom_vline(xintercept = 
                 mean(data[which(data$time == 0 & data[var_sex] == female_male[2]), var],
                      na.rm = TRUE),
               color = colors["male"], size = 1) +
    labs(
      title = paste("Distribution of", var_name, "<br>", color_text),
      x = var_name,
      y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_markdown()
    )
  return(p)
} # END plot_hist

# function for missingness analysis
miss_analysis <- function(ID, var, var_demo, var_rela, data, baseline = 0){
  missing <- data.frame(ID = data[ID])
  missing$missing <- apply(X = select(data, all_of(var)),
                           MARGIN = 1,
                           FUN = function(x) as.numeric(any(is.na(x))))
  missing <- missing %>% 
    group_by(!!sym(ID)) %>% 
    filter(missing == max(missing)) %>% 
    unique()
  
  # merge in demographic info
  demo <- data %>% filter(time == baseline) %>% select(all_of(c(ID, var_demo))) 
  missing <- merge(missing, demo)
  
  # merge in average scores
  rela <- data %>% filter(time == baseline) %>% select(all_of(c(ID, var_rela)))
  missing <- merge(missing, rela)      
  
  # return
  for(v in c(var_demo, var_rela)) {
    summary(
      glm(formula = paste("missing ~", v), family = "binomial", 
          data = missing))$coefficients %>%
      as.data.frame() %>%
      knitr::kable(
        caption = paste(v, "predicting missingness")) %>%
      kable_styling() %>%
      print()
  }
} # END miss_analysis

# function for test-retest reliability
icc_function <- function(ID, var, time_var, df){
  
  # select only variables of interest
  dat <- df %>%
    select(all_of(c(ID, var, time_var)))
  
  # pivot_wider
  dat <- dat %>%
    pivot_wider(names_from = time_var, values_from = var) %>%
    select(-all_of(c(ID)))
  
  # calculate ICC
  ICC_result = psych::ICC(x = dat, lmer = T)
  return(ICC_result$results %>% filter(type == "ICC1k"))
} # END ICC_function
