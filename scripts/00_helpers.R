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
  time_list <- unique(df$time)
  
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

h3_function <- function(var_list, df, baseline = 1) {
  
  # create cor_tab dataframe for correlated slopes results
  cor_tab <- data.frame()
  
  # create slopes_tab dataframe for fixed effects estimates
  slopes_tab <- data.frame()
  
  # create slope_df dataframe for extracted slopes from separate gender models
  slope_df <- data.frame(couple = unique(df$couple),
                         time = baseline)
  
  # loop through variables
  for (var in var_list) {
    
    # var name for each participant
    p1_var <- paste0(var, "_1")
    p2_var <- paste0(var, "_2")
    
    # fit linear mixed model
    m1 <- lmer(formula = paste(p1_var, "~ time + (1 + time | couple)"),
               control = lmerControl(optimizer ="Nelder_Mead"),
               data = df)
    m2 <- lmer(formula = paste(p2_var, "~ time + (1 + time | couple)"),
               control = lmerControl(optimizer ="Nelder_Mead"),
               data = df)
    
    # extract and store fixed effects estimates
    slopes_tab <- rbind(
      slopes_tab,
      c("female", var, round(summary(m1)$coefficients[2,], 3)),
      c("male", var, round(summary(m2)$coefficients[2,], 3))
    )
    
    # extract fitted slopes
    slope_1 <- coef(m1)$couple[,2]
    slope_2 <- coef(m2)$couple[,2]
    
    # store fitted slopes
    slope_df <- merge(slope_df,
                      data.frame(couple = rownames(coef(m1)$couple),
                                 slope_1,
                                 slope_2))
    names(slope_df)[names(slope_df) == "slope_1"] <- paste0("slope_",
                                                            var, "_1")
    names(slope_df)[names(slope_df) == "slope_2"] <- paste0("slope_",
                                                            var, "_2")
    
    # extract correlation and p-value
    current_cor <- cor.test(slope_1, slope_2, method = "pearson",
                            alternative = "two.sided")
    
    # store values
    cor_tab <- rbind(cor_tab, c(var, 
                                round(current_cor$estimate, 3), 
                                round(current_cor$p.value, 5),
                                round(current_cor$conf.int[1], 3),
                                round(current_cor$conf.int[2], 3)))  
  } # END for var LOOP
  
  # rename cor_tab and ensure variable type
  names(cor_tab) <- c("variable", "correlation", "p_value", "LL", "UL")
  cor_tab$correlation <- as.numeric(cor_tab$correlation)
  cor_tab$p_value <- as.numeric(cor_tab$p_value)
  cor_tab$LL <- as.numeric(cor_tab$LL)
  cor_tab$UL <- as.numeric(cor_tab$UL)
  
  # rename slopes_tab and ensure variable type
  names(slopes_tab) <- c("gender", "variable", 
                         "slope", "SE", "df", "t_value", "p_value")
  
  # return cor_tab and slope_df
  return(list(cor_tab = cor_tab,
              slopes_tab = slopes_tab,
              slope_df = slope_df))
} # END h3_function

# Research Question 2: Benefits of Assortative Mating --------------------------

## H4 Function: Baseline Benefits ----------------------------------------------

h4_function <- function(
    # character vector of personality variables
    var_list, 
    # character vector of profile names
    #   should match in df: 
    #     profilename_raw_r, profilename_centered_r, profilename_std_r
    prof_list,
    # character vector of relationship quality variables
    quality_list, 
    # time point defaults to baseline
    time = 1, df) {
  
  # create data frame to store results
  # INTERACTION_TAB for multiple regression
  #   nrow = # rela quality X # personality variables X 2 genders
  #   ncol = 12
  #     3 identifiers for quality and personality variables + gender
  #     3 variables for actor effect b, t, and p-value
  #     3 variables for partner effect b, t, and p-value
  #     3 variables for interaction effect b, t, and p-value
  interaction_tab <- data.frame()
  # DIFFERENCE_TAB for absolute difference simple regression
  #   nrow = # rela quality X # personality variables X 2 genders
  #   ncol = 6
  #     3 identifiers for quality and personality variables + gender
  #     3 variables for linear slope b, t, and p-value
  difference_tab <- data.frame()
  # PROFILE_TAB for multiple regression
  #   nrow = # rela quality X # personality profiles X 2 genders
  #   ncol = 6
  #     3 identifiers for quality and personality variables + gender
  #     9 variables for (linear slope b, t, and p-value) x 3 profiles raw, cent, std
  profile_tab <- data.frame()
  
  # iterate through relationship quality
  for (quality in quality_list) {
    
    # extract vector of each partner score on relationship quality variables
    p1_qual <- df[df$time == time, paste0(quality, "_1"), drop = T]
    p2_qual <- df[df$time == time, paste0(quality, "_2"), drop = T]
    
    # iterate through personality variables
    for (var in var_list) {

      # extract vector of each partner score on personality variables
      p1_var   <- df[df$time == time, paste0(var, "_1"), drop = T]
      p2_var   <- df[df$time == time, paste0(var, "_2"), drop = T]
      
      ### MULTIPLE REGRESSION ###
      
      # fit multiple regression models
      mod_1 <- summary(lm(p1_qual ~ p1_var * p2_var))$coefficients
      mod_2 <- summary(lm(p2_qual ~ p2_var * p1_var))$coefficients
      
      # extract coefficients
      actor_1   <- round(mod_1["p1_var", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      actor_2   <- round(mod_2["p2_var", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      partner_1 <- round(mod_1["p2_var", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      partner_2 <- round(mod_2["p1_var", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      int_1     <- round(mod_1["p1_var:p2_var", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      int_2     <- round(mod_2["p2_var:p1_var", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      
      # store values in interaction_tab
      interaction_tab <- rbind(
        interaction_tab,
        
        # first partner - female
        c(quality, var, "female", actor_1, partner_1, int_1),
        
        # second partner - male
        c(quality, var, "male", actor_2, partner_2, int_2)
      )
      
      ### ABSOLUTE DIFFERENCE ###
      
      # create a vector of absolute difference scores
      diff_var <- abs(p1_var - p2_var)
      
      # fit simple regression models
      mod_1 <- summary(lm(p1_qual ~ diff_var))$coefficients
      mod_2 <- summary(lm(p2_qual ~ diff_var))$coefficients
      
      # extract coefficients
      diff_1   <- round(mod_1["diff_var", 
                              c("Estimate", "t value", "Pr(>|t|)")], 3)
      diff_2   <- round(mod_2["diff_var", 
                              c("Estimate", "t value", "Pr(>|t|)")], 3)
      
      # store values in difference_tab
      difference_tab <- rbind(
        difference_tab,
        
        # first partner - female
        c(quality, var, "female", diff_1),
        
        # second partner - male
        c(quality, var, "male", diff_2)
      )
    } # END for var LOOP
    
    ### PROFILE REGRESSION ###
    
    # iterate through profile list
    for (prof in prof_list) {
      
      # extract profile correlations and transform to fisher's z scores
      profz_raw <- fisherz(
        rho = df[df$time == time, paste0(prof, "_raw_r"), drop = T]
      )
      profz_centered <- fisherz(
        rho = df[df$time == time, paste0(prof, "_centered_r"), drop = T]
      )
      profz_std <- fisherz(
        rho = df[df$time == time, paste0(prof, "_std_r"), drop = T]
      )
      
      # fix inifinity problem
      profz_raw[which(profz_raw == Inf)] <- 0.99999
      profz_centered[which(profz_centered == Inf)] <- 0.99999
      profz_std[which(profz_std == Inf)] <- 0.99999
      profz_raw[which(profz_raw == -Inf)] <- -0.99999
      profz_centered[which(profz_centered == -Inf)] <- -0.99999
      profz_std[which(profz_std == -Inf)] <- -0.99999
      
      # fit simple regression models
      mod_raw_1 <- summary(lm(p1_qual ~ profz_raw))$coefficients
      mod_raw_2 <- summary(lm(p2_qual ~ profz_raw))$coefficients
      mod_cen_1 <- summary(lm(p1_qual ~ profz_centered))$coefficients
      mod_cen_2 <- summary(lm(p2_qual ~ profz_centered))$coefficients
      mod_std_1 <- summary(lm(p1_qual ~ profz_std))$coefficients
      mod_std_2 <- summary(lm(p2_qual ~ profz_std))$coefficients
      
      # extract coefficients
      raw_1   <- round(mod_raw_1["profz_raw", 
                                c("Estimate", "t value", "Pr(>|t|)")], 3)
      raw_2   <- round(mod_raw_2["profz_raw", 
                                c("Estimate", "t value", "Pr(>|t|)")], 3)
      cen_1   <- round(mod_cen_1["profz_centered", 
                                c("Estimate", "t value", "Pr(>|t|)")], 3)
      cen_2   <- round(mod_cen_2["profz_centered", 
                                c("Estimate", "t value", "Pr(>|t|)")], 3)
      std_1   <- round(mod_std_1["profz_std", 
                                c("Estimate", "t value", "Pr(>|t|)")], 3)
      std_2   <- round(mod_std_2["profz_std", 
                                c("Estimate", "t value", "Pr(>|t|)")], 3)
      # store values in profile_tab
      profile_tab <- rbind(
        profile_tab,
        
        # first partner - female
        c(quality, prof, "female", raw_1, cen_1, std_1),
        
        # second partner - male
        c(quality, prof, "male", raw_2, cen_2, std_2)
      )
    }
    
  } # END for satis LOOP
  
  # rename interaction_tab columns
  names(interaction_tab) <- c("quality", "personality", "gender",
                              "actor_est", "actor_tval", "actor_pval",
                              "partner_est", "partner_tval", "partner_pval",
                              "int_est", "int_tval", "int_pval")
  # rename difference_tab columns
  names(difference_tab) <- c("quality", "personality", "gender",
                             "diff_est", "diff_tval", "diff_pval")
  
  # rename profile_tab columns
  names(profile_tab) <- c("quality", "profile", "gender",
                          "raw_est", "raw_tval", "raw_pval",
                          "cen_est", "cen_tval", "cen_pval",
                          "std_est", "std_tval", "std_pval")
  
  # return results
  return(list(interaction_tab = interaction_tab,
              difference_tab = difference_tab,
              profile_tab = profile_tab))
  
} # END h4_function

## H5 Function: Baseline Benefits with Longitudinal Predictors -----------------

h5_function <- function(
    # character vector of personality variables
    var_list, 
    # character vector of relationship quality variables
    quality_list, 
    # time point defaults to baseline
    time = 1, df) {
  
  # create data frame to store results
  # INTERACTION_TAB for multiple regression
  #   nrow = # rela quality X # personality variables X 2 genders
  #   ncol = 12
  #     3 identifiers for quality and personality variables + gender
  #     3 variables for actor effect b, t, and p-value
  #     3 variables for partner effect b, t, and p-value
  #     3 variables for interaction effect b, t, and p-value
  interaction_tab <- data.frame()
  # DIFFERENCE_TAB for absolute difference simple regression
  #   nrow = # rela quality X # personality variables X 2 genders
  #   ncol = 6
  #     3 identifiers for quality and personality variables + gender
  #     3 variables for linear slope b, t, and p-value
  difference_tab <- data.frame()

  # iterate through relationship quality
  for (quality in quality_list) {
    
    # extract vector of each partner score on relationship quality variables
    p1_qual <- df[df$time == time, paste0(quality, "_1"), drop = T]
    p2_qual <- df[df$time == time, paste0(quality, "_2"), drop = T]
    
    # iterate through personality variables
    for (var in var_list) {

      # extract vector of each partner slopes on personality variables
      p1_slope   <- df[df$time == time, paste0("slope_", var, "_1"), drop = T]
      p2_slope   <- df[df$time == time, paste0("slope_", var, "_2"), drop = T]
      
      ### MULTIPLE REGRESSION ###
      
      # fit multiple regression models
      mod_1 <- summary(lm(p1_qual ~ p1_slope * p2_slope))$coefficients
      mod_2 <- summary(lm(p2_qual ~ p2_slope * p1_slope))$coefficients
      
      # extract coefficients
      actor_1   <- round(mod_1["p1_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      actor_2   <- round(mod_2["p2_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      partner_1 <- round(mod_1["p2_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      partner_2 <- round(mod_2["p1_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      int_1     <- round(mod_1["p1_slope:p2_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      int_2     <- round(mod_2["p2_slope:p1_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      
      # store values in interaction_tab
      interaction_tab <- rbind(
        interaction_tab,
        
        # first partner - female
        c(quality, var, "female", actor_1, partner_1, int_1),
        
        # second partner - male
        c(quality, var, "male", actor_2, partner_2, int_2)
      )
      
      ### ABSOLUTE DIFFERENCE ###
      
      # create a vector of absolute difference scores
      diff_slope <- abs(p1_slope - p2_slope)
      
      # fit simple regression models
      mod_1 <- summary(lm(p1_qual ~ diff_slope))$coefficients
      mod_2 <- summary(lm(p2_qual ~ diff_slope))$coefficients
      
      # extract coefficients
      diff_1   <- round(mod_1["diff_slope", 
                              c("Estimate", "t value", "Pr(>|t|)")], 3)
      diff_2   <- round(mod_2["diff_slope", 
                              c("Estimate", "t value", "Pr(>|t|)")], 3)
      
      # store values in difference_tab
      difference_tab <- rbind(
        difference_tab,
        
        # first partner - female
        c(quality, var, "female", diff_1),
        
        # second partner - male
        c(quality, var, "male", diff_2)
      )
    } # END for var LOOP
    
  } # END for quality LOOP
    
  # rename interaction_tab columns
  names(interaction_tab) <- c("quality", "personality_slope", "gender",
                              "actor_est", "actor_tval", "actor_pval",
                              "partner_est", "partner_tval", "partner_pval",
                              "int_est", "int_tval", "int_pval")
  # rename difference_tab columns
  names(difference_tab) <- c("quality", "personality_slope", "gender",
                             "diff_est", "diff_tval", "diff_pval")
  
  # return results
  return(list(interaction_tab = interaction_tab,
              difference_tab = difference_tab))
} # END h5_function

## H6 Function: Longitudinal Benefits ------------------------------------------

h6_function <- function(
    # character vector of personality variables
    var_list, 
    # character vector of relationship quality variables
    quality_list, 
    # data frame for analyses
    df, baseline = 1) {
  
  # create data frame to store results
  # INTERACTION_TAB for multiple regression
  #   nrow = # rela quality X # personality variables X 2 genders
  #   ncol = 12
  #     3 identifiers for quality and personality variables + gender
  #     3 variables for actor effect b, t, and p-value
  #     3 variables for partner effect b, t, and p-value
  #     3 variables for interaction effect b, t, and p-value
  interaction_tab <- data.frame()
  # DIFFERENCE_TAB for absolute difference simple regression
  #   nrow = # rela quality X # personality variables X 2 genders
  #   ncol = 6
  #     3 identifiers for quality and personality variables + gender
  #     3 variables for linear slope b, t, and p-value
  difference_tab <- data.frame()

  # iterate through relationship quality
  for (quality in quality_list) {
    
    # var name for each participant relationship quality
    p1_qual <- paste0(quality, "_1")
    p2_qual <- paste0(quality, "_2")
    
    # fit linear mixed model
    m1 <- lmer(formula = paste(p1_qual, "~ time + (1 + time | couple)"),
               control = lmerControl(optimizer ="Nelder_Mead"),
               data = df)
    m2 <- lmer(formula = paste(p2_qual, "~ time + (1 + time | couple)"),
               control = lmerControl(optimizer ="Nelder_Mead"),
               data = df)
    
    # extract fitted slopes
    p1_slope_qual <- coef(m1)$couple[,2]
    p2_slope_qual <- coef(m2)$couple[,2]

    # iterate through personality variables
    for (var in var_list) {

      # extract vector of each partner slopes on personality variables
      p1_slope   <- df[df$time == baseline, paste0("slope_", var, "_1"), drop = T]
      p2_slope   <- df[df$time == baseline, paste0("slope_", var, "_2"), drop = T]
      
      ### MULTIPLE REGRESSION ###
      
      # fit multiple regression models
      mod_1 <- summary(lm(p1_slope_qual ~ p1_slope * p2_slope))$coefficients
      mod_2 <- summary(lm(p2_slope_qual ~ p2_slope * p1_slope))$coefficients
      
      # extract coefficients
      actor_1   <- round(mod_1["p1_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      actor_2   <- round(mod_2["p2_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      partner_1 <- round(mod_1["p2_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      partner_2 <- round(mod_2["p1_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      int_1     <- round(mod_1["p1_slope:p2_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      int_2     <- round(mod_2["p2_slope:p1_slope", 
                               c("Estimate", "t value", "Pr(>|t|)")], 3)
      
      # store values in interaction_tab
      interaction_tab <- rbind(
        interaction_tab,
        
        # first partner - female
        c(quality, var, "female", actor_1, partner_1, int_1),
        
        # second partner - male
        c(quality, var, "male", actor_2, partner_2, int_2)
      )
      
      ### ABSOLUTE DIFFERENCE ###
      
      # create a vector of absolute difference scores
      diff_slope <- abs(p1_slope - p2_slope)
      
      # fit simple regression models
      mod_1 <- summary(lm(p1_slope_qual ~ diff_slope))$coefficients
      mod_2 <- summary(lm(p2_slope_qual ~ diff_slope))$coefficients
      
      # extract coefficients
      diff_1   <- round(mod_1["diff_slope", 
                              c("Estimate", "t value", "Pr(>|t|)")], 3)
      diff_2   <- round(mod_2["diff_slope", 
                              c("Estimate", "t value", "Pr(>|t|)")], 3)
      
      # store values in difference_tab
      difference_tab <- rbind(
        difference_tab,
        
        # first partner - female
        c(quality, var, "female", diff_1),
        
        # second partner - male
        c(quality, var, "male", diff_2)
      )
    } # END for var LOOP
    
  } # END for quality LOOP
    
  # rename interaction_tab columns
  names(interaction_tab) <- c("quality_slope", "personality_slope", "gender",
                              "actor_est", "actor_tval", "actor_pval",
                              "partner_est", "partner_tval", "partner_pval",
                              "int_est", "int_tval", "int_pval")
  # rename difference_tab columns
  names(difference_tab) <- c("quality_slope", "personality_slope", "gender",
                             "diff_est", "diff_tval", "diff_pval")
  
  # return results
  return(list(interaction_tab = interaction_tab,
              difference_tab = difference_tab))
} # END h6_function

## H7 Function: Cross-Lagged Effects -------------------------------------------

h7_function <- function(
  # character vector of personality variables
  var_list, 
  # character vector of personality profiles
  prof_list,
  # character vector of relationship quality variables
  quality_list, 
  # data frame for analyses
  df) {
  
  # extract number of time points
  n_times <- length(unique(df$time))
  
  # initialize results df
  # est_df: standardized solutions for cross-lagged paths
  # estprof_df: standardized solutions for cross-lagged paths of profile z
  est_df <- estprof_df <- data.frame()
  # fit_df: fit statistics for each model
  # fitprof_df: fit statistics for each model of profile z
  fit_df <- fitprof_df <- data.frame()
  
  # duplicate dataframe for variable creation
  wide_dat <- df
  
  # create difference score at each time point for personality variables
  for(var in var_list){
    wide_dat <- wide_dat %>%
      mutate(!!sym(paste0(var, "_diff")) := 
               abs(!!sym(paste0(var, "_1")) - !!sym(paste0(var, "_2"))))
  }
  
  # create fisher z-transformed score of profile correlations
  for(prof in prof_list){
    wide_dat <- wide_dat %>%
      mutate(!!sym(paste0(prof, "_raw_z")) :=
               fisherz(rho = !!sym(paste0(prof, "_raw_r"))),
             !!sym(paste0(prof, "_centered_z")) :=
               fisherz(rho = !!sym(paste0(prof, "_centered_r"))),
             !!sym(paste0(prof, "_std_z")) :=
               fisherz(rho = !!sym(paste0(prof, "_std_r"))))
  }
  
  # create wide dataframe with a column for each timepoint
  wide_dat <- wide_dat %>%
    pivot_wider(id_cols = "couple",
                names_from = "time",
                names_sep = "_t",
                values_from = c(
                  # relationship quality
                  paste0(rep(quality_list, each = 2), "_", c(1:2)),
                  
                  # univariate similarity
                  paste0(var_list, "_diff"),
                  
                  # profile correlations fisher z-transformed
                  paste0(prof_list, "_", c("raw_z", "centered_z", "std_z")))
    )
  
  # create and store within-subject list for lavaan formula
  within_qual_list <- paste0("w_qual_t", 1:n_times)
  within_diff_list <- paste0("w_diff_t", 1:n_times)
  # analogous versions for profiles
  within_profile_list <- paste0("w_profz_t", 1:n_times)
  
  # create and store within-subject cross-lagged effects
  #   so a later time point is regressed on t-1 time point
  #   not vice-versa
  cross_lagged <- c()
  for(time in 2:n_times){
    cross_lagged <- append(cross_lagged,
                           c(
                             paste0("w_qual_t", time, " ~ w_qual_t", (time-1)),
                             paste0("w_diff_t", time, " ~ w_diff_t", (time-1)),
                             paste0("w_qual_t", time, " ~ w_diff_t", (time-1)),
                             paste0("w_diff_t", time, " ~ w_qual_t", (time-1))))
  } # END for time LOOP
  # analogous versions for profiles
  cross_lagged_profile <- gsub(x = cross_lagged, pattern = "_diff_",
                               replacement = "_profz_")
  
  # create and store within-subject covariances at a given time point
  within_covar <- c(paste0("w_qual_t", 1:n_times, " ~~ w_diff_t", 1:n_times))
  # analogous versions for profiles
  within_covar_profile <- gsub(x = within_covar, pattern = "_diff_",
                               replacement = "_profz_")
  
  # create and store within-subject variances at a given time point
  within_var <- c(paste0("w_qual_t", 1:n_times, " ~~ w_qual_t", 1:n_times),
                  paste0("w_diff_t", 1:n_times, " ~~ w_diff_t", 1:n_times))
  # analogous versions for profiles
  within_var_profile <- gsub(x = within_var, pattern = "_diff_",
                             replacement = "_profz_")
  
  # iterate through relationship quality
  for(qual in quality_list){
    
    ### UNIVARIATE SIMILARITY ###
    
    # iterate through personality variable
    for(var in var_list) {
      
      # store RI-CLPM model specifications
      #   model_p1: relationship quality reported by P1
      model_p1 <- paste(
      
        # random intercepts for relationship quality
        #   reported by P1
        #   RI_qual =~ 1*qual_1_t1 + 1*qual_1_t2 + ...
        "RI_qual =~ 1*", paste0(qual, "_1_t", 1:n_times, collapse = "+1*"),
        
        # random intercepts for personality similarity/difference
        #   RI_diff =~ 1*diff_t1 + 1*diff_t2 + ...
        "
        RI_diff =~ 1*", paste0(var, "_diff_t", 1:n_times, collapse = "+1*"),
        "
        ",
        
        # within-subject quality reported by P1 at each time point
        #   w_qual_t1 =~ 1*qual_1_t1
        #   w_qual_t2 =~ 1*qual_1_t2
        #   ...
        paste0(within_qual_list, " =~ 1*", qual, "_1_t", 1:n_times, 
               # each definition on 1 line
               collapse = " \n "),
        "
        ",
        
        # within-subject personality difference at each time point
        #   w_diff_t1 =~ 1*diff_t1
        #   w_diff_t2 =~ 1*diff_t2
        #   ...
        paste0(within_diff_list, " =~ 1*", var, "_diff_t", 1:n_times, 
               # each definition on 1 line
               collapse = " \n "),
        "
        ",
        
        # cross-lagged effects among within-subjects definitions
        # so a later time point is regressed on t-1 time point
        #   w_qual_t3 + w_diff_t3 ~ w_qual_t2 + w_diff_t2
        #   w_qual_t2 + w_diff_t2 ~ w_qual_t1 + w_diff_t1
        #   ...
        paste(cross_lagged, collapse = " \n "),
        "
        ",
        
        # covariances among within-subject effects
        #   w_qual_t1 ~~ w_diff_t1
        #   w_qual_t2 ~~ w_diff_t2
        #   ...
        paste(within_covar, collapse = " \n "),
        
        # covariances among random-intercepts
        "
        RI_qual ~~ RI_diff",
  
        #variances among random-intercepts
        "
        RI_qual ~~ RI_qual
        RI_diff ~~ RI_diff
        ",
        
        # residual variances among within-subject effects
        paste(within_var, collapse = " \n ")
      ) # END model_p1 definition
      
      # store RI-CLPM model specifications
      # model_p2: relationship quality reported by P2
      model_p2 <- gsub(x = model_p1, pattern = "_1_", replacement = "_2_")
      
      # fit lavaan model
      fit_p1 <- lavaan(model_p1, data = wide_dat, missing = "FIML", 
                       meanstructure = T, int.ov.free = T)
      fit_p2 <- lavaan(model_p2, data = wide_dat, missing = "FIML", 
                       meanstructure = T, int.ov.free = T)
      
      # extract fit statistics
      fitstats_p1 <- round(
        fitMeasures(fit_p1, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                    output = "vector"),
        3)
      fitstats_p2 <- round(
        fitMeasures(fit_p2, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                    output = "vector"),
        3)
      
      # store fit statistics
      fit_df <- rbind(fit_df,
                      c(1, var, qual, fitstats_p1),
                      c(2, var, qual, fitstats_p2))
      
      # standardized cross-lagged solutions
      est_p1 <- as.data.frame(unclass(standardizedSolution(fit_p1))) %>% 
        filter(op == "~") %>%
        mutate(partner = 1,
               personality = var,
               quality = qual) %>%
        relocate(partner, personality, quality) %>%
          mutate_if(is.numeric, function(x){round(x,3)})
      est_p2 <- as.data.frame(unclass(standardizedSolution(fit_p2))) %>% 
        filter(op == "~") %>%
        mutate(partner = 2,
               personality = var,
               quality = qual) %>%
        relocate(partner, personality, quality) %>%
          mutate_if(is.numeric, function(x){round(x,3)})
      
      # store solutions
      est_df <- rbind(est_df,
                      est_p1,
                      est_p2)
    } # END for var LOOP
    
    
    ### PROFILE SIMILARITY ###
    
    # iterate through personality profiles
    for(prof in prof_list) {
      
      # iterate through raw/centered/std profiles
      for(proftype in c("raw", "centered", "std")) {
        # store RI-CLPM model specifications
        #   modelprof_p1: relationship quality reported by P1
        modelprof_p1 <- paste(
        
          # random intercepts for relationship quality
          #   reported by P1
          #   RI_qual =~ 1*qual_1_t1 + 1*qual_1_t2 + ...
          "RI_qual =~ 1*", paste0(qual, "_1_t", 1:n_times, collapse = "+1*"),
          
          # random intercepts for profile correlation
          #   RI_profcor =~ 1*diff_t1 + 1*diff_t2 + ...
          "
          RI_profz =~ 1*", paste0(prof, "_", proftype, "_z_t", 1:n_times, 
                                  collapse = "+1*"),
          "
          ",
          
          # within-subject quality reported by P1 at each time point
          #   w_qual_t1 =~ 1*qual_1_t1
          #   w_qual_t2 =~ 1*qual_1_t2
          #   ...
          paste0(within_qual_list, " =~ 1*", qual, "_1_t", 1:n_times, 
                 # each definition on 1 line
                 collapse = " \n "),
          "
          ",
          
          # within-subject personality difference at each time point
          #   w_profz_t1 =~ 1*raw_z_t1
          #   w_profz_t2 =~ 1*raw_z_t2
          #   ...
          paste0(within_profile_list, " =~ 1*", 
                 prof, "_", proftype, "_z_t", 1:n_times, 
                 # each definition on 1 line
                 collapse = " \n "),
          "
          ",
          
          # cross-lagged effects among within-subjects definitions
          # so a later time point is regressed on t-1 time point
          #   w_qual_t3 + w_profz_t3 ~ w_qual_t2 + w_profz_t2
          #   w_qual_t2 + w_profz_t2 ~ w_qual_t1 + w_profz_t1
          #   ...
          paste(cross_lagged_profile, collapse = " \n "),
          "
          ",
          
          # covariances among within-subject effects
          #   w_qual_t1 ~~ w_profz_t1
          #   w_qual_t2 ~~ w_profz_t2
          #   ...
          paste(within_covar_profile, collapse = " \n "),
          
          # covariances among random-intercepts
          "
          RI_qual ~~ RI_profz",
    
          #variances among random-intercepts
          "
          RI_qual ~~ RI_qual
          RI_profz ~~ RI_profz
          ",
          
          # residual variances among within-subject effects
          paste(within_var_profile, collapse = " \n ")
        ) # END model_p1 definition
        
        # store RI-CLPM model specifications
        # modelprof_p2: relationship quality reported by P2
        modelprof_p2 <- gsub(x = modelprof_p1, pattern = "_1_", 
                             replacement = "_2_")
        
        # fit lavaan model
        fitprof_p1 <- lavaan(modelprof_p1, data = wide_dat, missing = "FIML", 
                         meanstructure = T, int.ov.free = T)
        fitprof_p2 <- lavaan(modelprof_p2, data = wide_dat, missing = "FIML", 
                         meanstructure = T, int.ov.free = T)
        
        # extract fit statistics
        fitstatsprof_p1 <- round(
          fitMeasures(fitprof_p1, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                      output = "vector"),
          3)
        fitstatsprof_p2 <- round(
          fitMeasures(fitprof_p2, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                      output = "vector"),
          3)
        
        # store fit statistics
        fitprof_df <- rbind(fitprof_df,
                            c(1, proftype, prof, qual, fitstatsprof_p1),
                            c(2, proftype, prof, qual, fitstatsprof_p2))
        
        # standardized cross-lagged solutions
        estprof_p1 <- as.data.frame(unclass(standardizedSolution(fitprof_p1))) %>% 
          filter(op == "~") %>%
          mutate(partner = 1,
                 type = proftype,
                 profile = prof,
                 quality = qual) %>%
          relocate(partner, type, profile, quality) %>%
          mutate_if(is.numeric, function(x){round(x,3)})

        estprof_p2 <- as.data.frame(unclass(standardizedSolution(fitprof_p2))) %>% 
          filter(op == "~") %>%
          mutate(partner = 2,
                 type = proftype,
                 profile = prof,
                 quality = qual) %>%
          relocate(partner, type, profile, quality) %>%
          mutate_if(is.numeric, function(x){round(x,3)})
        
        # store solutions
        estprof_df <- rbind(estprof_df,
                            estprof_p1,
                            estprof_p2)    
      } # END for proftype LOOP
    } # END for prof LOOP

  } # END for qual LOOP
  
  # returns results
  names(fit_df) <- c("partner", "personality", "quality",
                     "chisq", "df", "pvalue", "cfi", "rmsea", "srmr")
  names(fitprof_df) <- c("partner", "type", "profile", "quality",
                         "chisq", "df", "pvalue", "cfi", "rmsea", "srmr")
  return(list(est_df = est_df,
              fit_df = fit_df,
              estprof_df = estprof_df,
              fitprof_df = fitprof_df))
} # END h7_function

# Research Question 3: Actor/Partner/Perceived/Actual Similarity ---------------

## H8 Function: Actor/Partner Effect on Quality --------------------------------

h8_function <- function(
  # character vector of all personality variables  
  var_list,
  
  # character vector of all relationship quality variables
  quality_list,
  
  # time point (default to baseline) and analytic dataframe
  time = 1, df) {
  
  # create df to store results
  fit_df <- est_df <- data.frame()
  
  # loop through personality variables
  for(var in var_list){
    df <- df %>%
      # filter to only include the specified time point
      #   default to 1 for baseline
      filter(time == time) %>%

      # create difference score for personality variables
      mutate(!!sym(paste0(var, "_diff")) := 
               abs(!!sym(paste0(var, "_1")) - !!sym(paste0(var, "_2"))))
    
    # loop through relationship satisfaction variables
    for(qual in quality_list){
      
      model <- paste(
        # female actor effect
        paste0(qual, "_1", " ~ a1*", var, "_1"),
        
        # male actor effect
        paste0(qual, "_2", " ~ a2*", var, "_2"),
        
        # female to male partner effect
        paste0(qual, "_2", " ~ p12*", var, "_1"),
        
        # male to female partner effect
        paste0(qual, "_1", " ~ p21*", var, "_2"),

        # dissimilarity to female effect
        paste0(qual, "_1", " ~ d1*", var, "_diff"),
        
        # dissimilarity to male effect
        paste0(qual, "_2", " ~ d2*", var, "_diff"),
        
        # female mean personality
        paste0(var, "_1", " ~ m1*1"),
        
        # male mean personality
        paste0(var, "_2", " ~ m2*1"),
        
        # dissimilarity mean
        paste0(var, "_diff", " ~ md*1"),
        
        # female intercept
        paste0(qual, "_1", " ~ i1*1"),
        
        # male intercept
        paste0(qual, "_2", " ~ i2*1"),
          
        # female variance for personality
        paste0(var, "_1", " ~~ v1*", var, "_1"),
        
        # male variance for personality
        paste0(var, "_2", " ~~ v2*", var, "_2"),

        # dissimilarity variance
        paste0(var, "_diff", " ~~ vd*", var, "_diff"),
        
        # female error variance for relationship quality
        paste0(qual, "_1", " ~~ e1*", qual, "_1"),

        # male error variance for relationship quality
        paste0(qual, "_2", " ~~ e2*", qual, "_2"),
        
        # covariance between female and male personality
        paste0(var, "_1", " ~~ v12*", var, "_2"),
          
        # covariance between female personality and dissimilarity
        paste0(var, "_1", " ~~ v1d*", var, "_diff"),
        
        # covariance between male personality and dissimilarity
        paste0(var, "_2", " ~~ v2d*", var, "_diff"),
        
        # error covariance between female and male relationship quality
        paste0(qual, "_1", " ~~ e12*", qual, "_2"),
        
        # each line is separated by a line break
        sep = " \n "
      ) # END model DEF
      
      # fit lavaan model
      fit <- lavaan(model, data = df, missing = "FIML")
      
      # standardized actor/partner/similarity solutions
      est <- as.data.frame(unclass(standardizedSolution(fit))) %>% 
        filter(op == "~") %>%
         mutate_if(is.numeric, function(x){round(x,3)})
      
      # store results in results df
      est_df <- rbind(
        est_df, est
      )

    } # END for qual LOOP
  } # END for var LOOP
  
  # return results
  return(est_df)
} # END h8_function

## H9 Function: Perceived/Actual Similarity Comparison -------------------------

h9_function <- function(
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
} # END h9_function

## H10 Function: Perceived/Actual Similarity Effect on Quality -----------------

h10_function <- function(
  # character vector of all personality variables with self/other reports
  #   naming convention: var_self and var_partner
  perception_list,
  
  # character vector of all relationship quality variables
  quality_list,
  
  # time point (default to baseline) and analytic dataframe
  time = 1, df) {

  # create dataframes to store results
  #   est_df: standardized solutions from APIM models
  #           separate for p1 and p2's perceptions
  est_df_p1 <- est_df_p2 <- data.frame()
  
  # loop through self/other variables
  for (var in perception_list) {
    
    # create difference scores for personality variables
    df <- df %>%
      # filter to only include the specified time point
      #   default to 1 for baseline
      filter(time == time) %>%

      # var_diff_p1: perceived dissimilarity between p1's self-perception 
      #              and p1's perception of p2
      # var_diff_p2: perceived dissimilarity between p2's self-perception
      #              and p2's perception of p1
      mutate(!!sym(paste0(var, "_diff_p1")) :=
               abs(!!sym(paste0(var, "_self_1")) - !!sym(paste0(var, "_partner_1"))),
             !!sym(paste0(var, "_diff_p2")) :=
               abs(!!sym(paste0(var, "_self_2")) - !!sym(paste0(var, "_partner_2"))))
    
    # loop through relationship satisfaction variables
    for(qual in quality_list){
      
      # model_p1 includes perceptions of female partner p1
      model_p1 <- paste(
        # female actor effect
        paste0(qual, "_1", " ~ a1*", var, "_self_1"),
        
        # female-perception actor effect
        paste0(qual, "_2", " ~ ap1*", var, "_partner_1"),
        
        # female to male partner effect
        paste0(qual, "_2", " ~ p12*", var, "_self_1"),
        
        # female-perception to female partner effect
        paste0(qual, "_1", " ~ p21*", var, "_partner_1"),

        # dissimilarity to female effect
        paste0(qual, "_1", " ~ d1*", var, "_diff_p1"),
        
        # dissimilarity to male effect
        paste0(qual, "_2", " ~ d2*", var, "_diff_p1"),
        
        # female mean personality
        paste0(var, "_self_1", " ~ m1*1"),
        
        # female-perception mean personality
        paste0(var, "_partner_1", " ~ mp1*1"),
        
        # dissimilarity mean
        paste0(var, "_diff_p1", " ~ md1*1"),
        
        # female intercept
        paste0(qual, "_1", " ~ i1*1"),
        
        # male intercept
        paste0(qual, "_2", " ~ i2*1"),
          
        # female variance for personality
        paste0(var, "_self_1", " ~~ v1*", var, "_self_1"),
        
        # female-perception variance for personality
        paste0(var, "_partner_1", " ~~ vp1*", var, "_partner_1"),

        # dissimilarity variance
        paste0(var, "_diff_p1", " ~~ vd1*", var, "_diff_p1"),
        
        # female error variance for relationship quality
        paste0(qual, "_1", " ~~ e1*", qual, "_1"),

        # male error variance for relationship quality
        paste0(qual, "_2", " ~~ e2*", qual, "_2"),
        
        # covariance between female and female-perception personality
        paste0(var, "_self_1", " ~~ v1p1*", var, "_partner_1"),
          
        # covariance between female personality and dissimilarity
        paste0(var, "_self_1", " ~~ v1d1*", var, "_diff_p1"),
        
        # covariance between female-perception personality and dissimilarity
        paste0(var, "_partner_1", " ~~ vp1d1*", var, "_diff_p1"),
        
        # error covariance between female and male relationship quality
        paste0(qual, "_1", " ~~ e12*", qual, "_2"),
        
        # each line is separated by a line break
        sep = " \n "
      ) # END model_p1 DEF
      
      # model_p1 includes perceptions of female partner p1
      model_p2 <- paste(
        
        # male actor effect
        paste0(qual, "_2", " ~ a2*", var, "_self_2"),
        
        # male-perception actor effect
        paste0(qual, "_1", " ~ ap2*", var, "_partner_2"),
        
        # male to female partner effect
        paste0(qual, "_1", " ~ p21*", var, "_self_2"),
        
        # male-perception to male partner effect
        paste0(qual, "_2", " ~ p21*", var, "_partner_2"),

        # dissimilarity to female effect
        paste0(qual, "_1", " ~ d1*", var, "_diff_p2"),
        
        # dissimilarity to male effect
        paste0(qual, "_2", " ~ d2*", var, "_diff_p2"),
        
        # male mean personality
        paste0(var, "_self_2", " ~ m2*1"),
        
        # male-perception mean personality
        paste0(var, "_partner_2", " ~ mp2*1"),
        
        # dissimilarity mean
        paste0(var, "_diff_p2", " ~ md2*1"),
        
        # female intercept
        paste0(qual, "_1", " ~ i1*1"),
        
        # male intercept
        paste0(qual, "_2", " ~ i2*1"),
          
        # male variance for personality
        paste0(var, "_self_2", " ~~ v2*", var, "_self_2"),
        
        # male-perception variance for personality
        paste0(var, "_partner_2", " ~~ vp2*", var, "_partner_2"),

        # dissimilarity variance
        paste0(var, "_diff_p2", " ~~ vd2*", var, "_diff_p2"),
        
        # female error variance for relationship quality
        paste0(qual, "_1", " ~~ e1*", qual, "_1"),

        # male error variance for relationship quality
        paste0(qual, "_2", " ~~ e2*", qual, "_2"),
        
        # covariance between male and male-perception personality
        paste0(var, "_self_2", " ~~ v2p2*", var, "_partner_2"),
          
        # covariance between male personality and dissimilarity
        paste0(var, "_self_2", " ~~ v2d2*", var, "_diff_p2"),
        
        # covariance between male-perception personality and dissimilarity
        paste0(var, "_partner_2", " ~~ vp2d2*", var, "_diff_p2"),
        
        # error covariance between female and male relationship quality
        paste0(qual, "_1", " ~~ e12*", qual, "_2"),
        
        # each line is separated by a line break
        sep = " \n "
      ) # END model_p2 DEF
      
      # fit lavaan model
      fit_p1 <- lavaan(model_p1, data = df, missing = "FIML")
      fit_p2 <- lavaan(model_p2, data = df, missing = "FIML")

      # standardized actor/partner/similarity solutions
      est_p1 <- as.data.frame(unclass(standardizedSolution(fit_p1))) %>% 
        filter(op == "~") %>%
        mutate_if(is.numeric, function(x){round(x,3)}) %>%
        mutate("perception" = "female") %>%
        relocate(perception)
      est_p2 <- as.data.frame(unclass(standardizedSolution(fit_p2))) %>% 
        filter(op == "~") %>%
        mutate_if(is.numeric, function(x){round(x,3)}) %>%
        mutate("perception" = "male") %>%
        relocate(perception)
      
      # store results
      est_df_p1 <- rbind(
        est_df_p1, est_p1
      )
      est_df_p2 <- rbind(
        est_df_p2, est_p2
      )
    } # END for qual LOOP
  } # END for var LOOP
  
  # return results
  return(list(est_df_p1 = est_df_p1,
              est_df_p2 = est_df_p2))
} # END h10_function

# Miscellaneous Functions ------------------------------------------------------


