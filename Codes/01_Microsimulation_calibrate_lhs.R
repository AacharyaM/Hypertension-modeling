#######################################################################
# Filename: 01_Microsimulation_calibrate_lhs.R
# Created by: Mahip
# Date created: 04/08/2025
# Purpose: Creating a Markov microsimulation model for
# patients >=60 years and calibrating the parameters
# Using Latin hypercube sampling with 1,000 samples for calibrating the model
# based on incidence of MI, HF, Stroke, and TIA; overall mortality, cv mortality
# and non-cv mortality, stroke mortality, and MI mortality
#######################################################################

# Loading all the required packages (installing them if they are not installed already)
loadpackages <- function(packagename) {
  if (!suppressWarnings(require(packagename, character.only = TRUE))) {
    
    # Install the package if not installed already
    install.packages(packagename)
    
    # Load the package
    library(packagename, character.only = TRUE)
  } else {
    library(packagename, character.only = TRUE)
  }
}

# Loading the individual libraries
loadpackages("here")
loadpackages("lhs")
loadpackages("parallel")
loadpackages("foreach")
loadpackages("doRNG")

# Also installing readxl, truncnorm and doParallel packages but not calling the entire library and only using where needed
installonly <- function(packagename) {
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
}

# Installing individual packages
installonly("readxl")
installonly("truncnorm")
installonly("doParallel")

# Changing printing formats
options(max.print = 1000000)
options(scipen = 999)

# Setting seed
set.seed(999)

# All parameters####
l_hrparams <- list(
  hr_mi = 3.19,
  hr_hf = 3.87,
  hr_stroke = 2.30,
  hr_tia = 1.35,
  hr_othercvd = 3.72,
  hr_chroniccvd = 1.35,
  hr_renal = 1.20,
  hr_advrenal = 4.03
)

# Number of cycles
n_cycles <- 612

# Number of states
n_states <- 71

# Vector of state names
cond1 <- c(
  "HTN", "MI", "MI2", "MI3", "MI4", "MI5", "MI6", "MI7", "MI8",
  "MI9", "MI10", "MI11", "MI12", "HF", "HF2", "HF3", "HF4", "HF5",
  "HF6", "HF7", "HF8", "HF9", "HF10", "HF11", "HF12", "Stroke", "Stroke2",
  "Stroke3", "Stroke4", "Stroke5", "Stroke6", "Stroke7", "Stroke8", "Stroke9",
  "Stroke10", "Stroke11", "Stroke12", "TIA", "TIA2", "TIA3", "TIA4", "TIA5",
  "TIA6", "TIA7", "TIA8", "TIA9", "TIA10", "TIA11", "TIA12", "Othercvd", "Othercvd2",
  "Othercvd3", "Othercvd4", "Othercvd5", "Othercvd6", "Othercvd7", "Othercvd8", "Othercvd9",
  "Othercvd10", "Othercvd11", "Othercvd12",
  "Renal", "Chroniccvd", "Advrenal", "MIdeath", "HFdeath", "Strokedeath", "Othercvddeath", "Noncvdeath", "Chroniccvddeath", "Death"
)

# Creating an array for transition probabilities across cycles

a_transition_matrix <- array(NA_real_, dim = c(n_states, n_states, n_cycles), dimnames = list(
  from = cond1,
  to = cond1,
  cycle = 1:n_cycles
))

# Importing an Excel file with annual transition rates
transition_rates <- readxl::read_excel("Data/Transition prob_for R_2021.xlsx", sheet = "rates")

# Converting the annual rates to monthly rates and then to monthly probabilities
transition_rates$mnth_rate <- (transition_rates$annual_rate) / 12
transition_rates$mnth_prob <- 1 - exp(-(transition_rates$annual_rate) / 12)

# Creating an indicator for the type of transition probability
transition_rates$prob_type <- paste("p", transition_rates$From, transition_rates$To, sep = "_")

# Importing background mortality info (source: CDC Life Table)
background_mortality <- readxl::read_excel("Data/Background mortality valid_for R.xlsx", sheet = "death")
background_mortality$mnth_rate <- (-log(1 - background_mortality$annual_prob)) / 12
background_mortality$mnth_prob <- 1 - exp(-(background_mortality$mnth_rate))

# Range of input parameters for changing in calibration
lowbound <- c(
  back_multiplier = 0.33,
  p_HTN_MI = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_MI")],
  p_HTN_HF = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_HF")],
  p_HTN_Stroke = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Stroke")],
  p_HTN_TIA = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_TIA")],
  p_HTN_Othercvd = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Othercvd")],
  p_MI_MI = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_MI")],
  p_MI_HF = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_HF")],
  p_MI_Stroke = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Stroke")],
  p_MI_TIA = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_TIA")],
  p_MI_Othercvd = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Othercvd")],
  p_HF_MI = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_MI")],
  p_HF_HF = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_HF")],
  p_HF_Stroke = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Stroke")],
  p_HF_TIA = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_TIA")],
  p_HF_Othercvd = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Othercvd")],
  p_Stroke_MI = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_MI")],
  p_Stroke_HF = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_HF")],
  p_Stroke_Stroke = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Stroke")],
  p_Stroke_TIA = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_TIA")],
  p_Stroke_Othercvd = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Othercvd")],
  p_TIA_MI = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_MI")],
  p_TIA_HF = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_HF")],
  p_TIA_Stroke = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Stroke")],
  p_TIA_TIA = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_TIA")],
  p_TIA_Othercvd = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Othercvd")],
  p_Othercvd_MI = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_MI")],
  p_Othercvd_HF = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_HF")],
  p_Othercvd_Stroke = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Stroke")],
  p_Othercvd_TIA = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_TIA")],
  p_Othercvd_Othercvd = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Othercvd")],
  p_Chroniccvd_MI = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_MI")],
  p_Chroniccvd_HF = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_HF")],
  p_Chroniccvd_Stroke = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_Stroke")],
  p_Chronicvd_TIA = 0.8 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_TIA")],
  hr_multiplier = 1.0
)

upperbound <- c(
  back_multiplier = 0.75,
  p_HTN_MI = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_MI")],
  p_HTN_HF = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_HF")],
  p_HTN_Stroke = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Stroke")],
  p_HTN_TIA = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_TIA")],
  p_HTN_Othercvd = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Othercvd")],
  p_MI_MI = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_MI")],
  p_MI_HF = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_HF")],
  p_MI_Stroke = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Stroke")],
  p_MI_TIA = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_TIA")],
  p_MI_Othercvd = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Othercvd")],
  p_HF_MI = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_MI")],
  p_HF_HF = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_HF")],
  p_HF_Stroke = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Stroke")],
  p_HF_TIA = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_TIA")],
  p_HF_Othercvd = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Othercvd")],
  p_Stroke_MI = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_MI")],
  p_Stroke_HF = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_HF")],
  p_Stroke_Stroke = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Stroke")],
  p_Stroke_TIA = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_TIA")],
  p_Stroke_Othercvd = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Othercvd")],
  p_TIA_MI = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_MI")],
  p_TIA_HF = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_HF")],
  p_TIA_Stroke = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Stroke")],
  p_TIA_TIA = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_TIA")],
  p_TIA_Othercvd = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Othercvd")],
  p_Othercvd_MI = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_MI")],
  p_Othercvd_HF = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_HF")],
  p_Othercvd_Stroke = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Stroke")],
  p_Othercvd_TIA = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_TIA")],
  p_Othercvd_Othercvd = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Othercvd")],
  p_Chroniccvd_MI = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_MI")],
  p_Chroniccvd_HF = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_HF")],
  p_Chroniccvd_Stroke = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_Stroke")],
  p_Chronicvd_TIA = 1.2 * transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_TIA")],
  hr_multiplier = 1.4
)

# Sampling the values for the 36 parameters
v_lhs_param_names <- c(
  "back_multiplier", "p_HTN_MI", "p_HTN_HF", "p_HTN_Stroke", "p_HTN_TIA", "p_HTN_Othercvd",
  "p_MI_MI", "p_MI_HF", "p_MI_Stroke", "p_MI_TIA", "p_MI_Othercvd",
  "p_HF_MI", "p_HF_HF", "p_HF_Stroke", "p_HF_TIA", "p_HF_Othercvd",
  "p_Stroke_MI", "p_Stroke_HF", "p_Stroke_Stroke", "p_Stroke_TIA", "p_Stroke_Othercvd",
  "p_TIA_MI", "p_TIA_HF", "p_TIA_Stroke", "p_TIA_TIA", "p_TIA_Othercvd",
  "p_Othercvd_MI", "p_Othercvd_HF", "p_Othercvd_Stroke", "p_Othercvd_TIA", "p_Othercvd_Othercvd",
  "p_Chroniccvd_MI", "p_Chroniccvd_HF", "p_Chroniccvd_Stroke", "p_Chroniccvd_TIA", "hr_multiplier"
)
n_lhs_params <- length(v_lhs_param_names)
n_lhs_samples <- 1000
m_lhs_space <- randomLHS(n_lhs_samples, n_lhs_params)

# Rescaling the parameters using lhs values
m_lhs_param_values <- matrix(nrow = n_lhs_samples, ncol = n_lhs_params)
for (i in 1:n_lhs_params) {
  m_lhs_param_values[, i] <- qunif(m_lhs_space[, i],
    min = lowbound[i],
    max = upperbound[i]
  )
}
# Adding the column names
colnames(m_lhs_param_values) <- v_lhs_param_names

# Calibration targets
n_lhs_targets <- 8
v_lhs_targets <- c("n_Death", "n_MI", "n_HF", "n_Stroke", "n_TIA", "n_Cvdeath", "n_Noncvdeath", "n_Strokedeath")
m_gof <- matrix(nrow = n_lhs_samples, ncol = n_lhs_targets * 2)

colnames(m_gof)[9:16] <- paste0(v_lhs_targets, "_fit")

colnames(m_gof)[1:8] <- v_lhs_targets

# Values for calibration targets
m_gof[, "n_Death"] <- 45000 # 169 cycles
m_gof[, "n_Cvdeath"] <- 24300 # 168 cycles
m_gof[, "n_Noncvdeath"] <- 20700 # 168 cycles
m_gof[, "n_Strokedeath"] <- 5000 # 264 cycles
m_gof[, "n_MI"] <- 5054 # 60 cycles #calculated using total number of stroke (159) and MI (98) and stroke rate of 8.2 per 100 in placebo
m_gof[, "n_HF"] <- 5483 # 60 cycles #calculated using total number of stroke reported in initial report (149 + 14), number of left ventricular failure in the same report (102 + 7) and stroke rate
m_gof[, "n_Stroke"] <- 8200 # 60 cycles
m_gof[, "n_TIA"] <- 4229 # 60 cycles  #calculated using total number of stroke (159) and TIA (82) events and stroke rate of 8.2 per 100 in placebo

# Creating a vector of numbers corresponding to disease states
cond1.index <- c(1:71)

# Number of individuals
n.i <- 100000

# Creating an empty array for the simulation cohort
a_state_membership <- array(NA_real_,
  dim = c(n_cycles, n_states),
  dimnames = list(
    cycle = 1:n_cycles,
    state = cond1
  )
)

# Creating a vector of n.i individuals with age drawn from a truncated
# normal distribution
v.age <- truncnorm::rtruncnorm(n.i, a = 60, b = 90, mean = 71.5, sd = 6.7)

# Age group 60 - 69 years : 41.8%
# Age group 70 - 79 years: 44.7%
# Age group 80 years or above: 13.5%

prop.age60_69 <- 0.418
prop.age70_79 <- 0.447
prop.age80_above <- 0.135
v.age60_69 <- v.age[v.age >= 60 & v.age <= 69.99]
v.age70_79 <- v.age[v.age >= 70 & v.age <= 79.99]
v.age80_above <- v.age[v.age >= 80 & v.age <= 90]

v.age60_69a <- sample(v.age60_69, size = round(prop.age60_69 * length(v.age)), replace = TRUE)
v.age70_79a <- sample(v.age70_79, size = round(prop.age70_79 * length(v.age)), replace = TRUE)
v.age80_abovea <- sample(v.age80_above, size = round(prop.age80_above * length(v.age)), replace = TRUE)

# Combining the age groups
v.age1 <- c(v.age60_69a, v.age70_79a, v.age80_abovea)

# Creating a vector with number of months from 60 years of age
v.age.monthsfrom60 <- floor((v.age1 - 60) * 12) + 1

# Creating a matrix for states
m.prob <- matrix(
  nrow = n.i, ncol = n_cycles,
  dimnames = list(
    paste("ind", 1:n.i, sep = " "),
    paste("cycle", 1:n_cycles, sep = " ")
  )
)

# Initializing the first cycle as HTN state
m.prob[, 1] <- "HTN"

# Using parallel processing to run the iterations
numclusters <- makeCluster(10)

# Moving the data frames and objects to all clusters
clusterExport(numclusters, c("m_lhs_param_values", "transition_rates", "l_hrparams", "a_state_membership", "background_mortality", "m_gof", "m.prob", "v.age.monthsfrom60"))

# Register the clusters
doParallel::registerDoParallel(numclusters)

# Using the random search approach for calibration
df_lhscalibration_results <- foreach(cb = 1:1000, .combine = rbind) %dorng% {
  # Function for transition probability####
  prob_func <- function(params) {
    with(params, {
      # Replacing the default values by values from each sample

      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_MI")] <- m_lhs_param_values[cb, "p_HTN_MI"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_HF")] <- m_lhs_param_values[cb, "p_HTN_HF"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Stroke")] <- m_lhs_param_values[cb, "p_HTN_Stroke"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_TIA")] <- m_lhs_param_values[cb, "p_HTN_TIA"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Othercvd")] <- m_lhs_param_values[cb, "p_HTN_Othercvd"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_MI")] <- m_lhs_param_values[cb, "p_MI_MI"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_HF")] <- m_lhs_param_values[cb, "p_MI_HF"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Stroke")] <- m_lhs_param_values[cb, "p_MI_Stroke"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_TIA")] <- m_lhs_param_values[cb, "p_MI_TIA"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Othercvd")] <- m_lhs_param_values[cb, "p_MI_Othercvd"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_MI")] <- m_lhs_param_values[cb, "p_HF_MI"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_HF")] <- m_lhs_param_values[cb, "p_HF_HF"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Stroke")] <- m_lhs_param_values[cb, "p_HF_Stroke"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_TIA")] <- m_lhs_param_values[cb, "p_HF_TIA"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Othercvd")] <- m_lhs_param_values[cb, "p_HF_Othercvd"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_MI")] <- m_lhs_param_values[cb, "p_Stroke_MI"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_HF")] <- m_lhs_param_values[cb, "p_Stroke_HF"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Stroke")] <- m_lhs_param_values[cb, "p_Stroke_Stroke"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_TIA")] <- m_lhs_param_values[cb, "p_Stroke_TIA"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Othercvd")] <- m_lhs_param_values[cb, "p_Stroke_Othercvd"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_MI")] <- m_lhs_param_values[cb, "p_TIA_MI"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_HF")] <- m_lhs_param_values[cb, "p_TIA_HF"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Stroke")] <- m_lhs_param_values[cb, "p_TIA_Stroke"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_TIA")] <- m_lhs_param_values[cb, "p_TIA_TIA"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Othercvd")] <- m_lhs_param_values[cb, "p_TIA_Othercvd"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_MI")] <- m_lhs_param_values[cb, "p_Othercvd_MI"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_HF")] <- m_lhs_param_values[cb, "p_Othercvd_HF"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Stroke")] <- m_lhs_param_values[cb, "p_Othercvd_Stroke"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_TIA")] <- m_lhs_param_values[cb, "p_Othercvd_TIA"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Othercvd")] <- m_lhs_param_values[cb, "p_Othercvd_Othercvd"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_MI")] <- m_lhs_param_values[cb, "p_Chroniccvd_MI"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_HF")] <- m_lhs_param_values[cb, "p_Chroniccvd_HF"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_Stroke")] <- m_lhs_param_values[cb, "p_Chroniccvd_Stroke"]
      transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_TIA")] <- m_lhs_param_values[cb, "p_Chroniccvd_TIA"]


      hr_mi <- hr_mi * m_lhs_param_values[cb, "hr_multiplier"]
      hr_hf <- hr_hf * m_lhs_param_values[cb, "hr_multiplier"]
      hr_stroke <- hr_stroke * m_lhs_param_values[cb, "hr_multiplier"]
      hr_tia <- hr_tia * m_lhs_param_values[cb, "hr_multiplier"]
      hr_othercvd <- hr_othercvd * m_lhs_param_values[cb, "hr_multiplier"]
      hr_chroniccvd <- hr_chroniccvd * m_lhs_param_values[cb, "hr_multiplier"]

      ## Filling the transition probability values####

      # From HTN to Death across all cycles
      # HTN to non-cardiovascular death
      for (x in 1:n_cycles) {
        a_transition_matrix["HTN", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HTN", "Noncvdeath", x] > 1] <- 1
      }

      # From HTN to Other states
      a_transition_matrix["HTN", "MI", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_MI")]
      a_transition_matrix["HTN", "HF", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_HF")]
      a_transition_matrix["HTN", "Stroke", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Stroke")]
      a_transition_matrix["HTN", "TIA", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_TIA")]
      a_transition_matrix["HTN", "Renal", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Renal")]
      a_transition_matrix["HTN", "Othercvd", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Othercvd")]
      a_transition_matrix["HTN", "HTN", ] <- 1 - apply(a_transition_matrix["HTN", , ], 2, sum, na.rm = TRUE)

      # From MI to all other states
      # MI to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI", "Noncvdeath", x] > 1] <- 1
      }

      # MI to MI2
      a_transition_matrix["MI", "MI2", ] <- 1 - apply(a_transition_matrix["MI", , ], 2, sum, na.rm = TRUE)

      # From MI2 to Other states
      # MI2 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI2", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI2 to Non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI2", "Noncvdeath", x] > 1] <- 1
      }

      # MI2 to MI3
      a_transition_matrix["MI2", "MI3", ] <- 1 - apply(a_transition_matrix["MI2", , ], 2, sum, na.rm = TRUE)

      # From MI3 to Other states
      # MI3 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI3", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI3", "Noncvdeath", x] > 1] <- 1
      }

      # MI3 to MI4
      a_transition_matrix["MI3", "MI4", ] <- 1 - apply(a_transition_matrix["MI3", , ], 2, sum, na.rm = TRUE)

      # From MI4 to Other states
      # MI4 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI4", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI4", "Noncvdeath", x] > 1] <- 1
      }

      # MI4 to MI5
      a_transition_matrix["MI4", "MI5", ] <- 1 - apply(a_transition_matrix["MI4", , ], 2, sum, na.rm = TRUE)

      # From MI5 to Other states
      # MI5 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI5", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI5", "Noncvdeath", x] > 1] <- 1
      }

      # MI5 to MI6
      a_transition_matrix["MI5", "MI6", ] <- 1 - apply(a_transition_matrix["MI5", , ], 2, sum, na.rm = TRUE)

      # From MI6 to Other states
      # MI6 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI6", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI6", "Noncvdeath", x] > 1] <- 1
      }

      # MI6 to MI7
      a_transition_matrix["MI6", "MI7", ] <- 1 - apply(a_transition_matrix["MI6", , ], 2, sum, na.rm = TRUE)

      # From MI7 to Other states
      # MI7 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI7", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI7", "Noncvdeath", x] > 1] <- 1
      }

      # MI7 to MI8
      a_transition_matrix["MI7", "MI8", ] <- 1 - apply(a_transition_matrix["MI7", , ], 2, sum, na.rm = TRUE)

      # From MI8 to Other states
      # MI8 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI8", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI8", "Noncvdeath", x] > 1] <- 1
      }

      # MI8 to MI9
      a_transition_matrix["MI8", "MI9", ] <- 1 - apply(a_transition_matrix["MI8", , ], 2, sum, na.rm = TRUE)

      # From MI9 to Other states
      # MI9 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI9", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI9", "Noncvdeath", x] > 1] <- 1
      }

      # MI9 to MI10
      a_transition_matrix["MI9", "MI10", ] <- 1 - apply(a_transition_matrix["MI9", , ], 2, sum, na.rm = TRUE)

      # From MI10 to Other states
      # MI10 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI10", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI10", "Noncvdeath", x] > 1] <- 1
      }

      # MI10 to MI11
      a_transition_matrix["MI10", "MI11", ] <- 1 - apply(a_transition_matrix["MI10", , ], 2, sum, na.rm = TRUE)

      # From MI11 to Other states
      # MI11 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI11", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI11", "Noncvdeath", x] > 1] <- 1
      }

      # MI11 to MI12
      a_transition_matrix["MI11", "MI12", ] <- 1 - apply(a_transition_matrix["MI11", , ], 2, sum, na.rm = TRUE)

      # From MI12 to Other states
      # MI12 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI12", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["MI12", "MIdeath", x] > 1] <- 1
      }

      # MI12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["MI12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["MI12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix["MI12", "MI", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_MI")]
      a_transition_matrix["MI12", "HF", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_HF")]
      a_transition_matrix["MI12", "Stroke", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Stroke")]
      a_transition_matrix["MI12", "TIA", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_TIA")]
      a_transition_matrix["MI12", "Othercvd", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_MI_Othercvd")]
      a_transition_matrix["MI12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix["MI12", , ], 2, sum, na.rm = TRUE)

      # From HF to Other states
      # HF to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF", "HFdeath", x] > 1] <- 1
      }

      # HF to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF", "Noncvdeath", x] > 1] <- 1
      }

      # HF to HF2
      a_transition_matrix["HF", "HF2", ] <- 1 - apply(a_transition_matrix["HF", , ], 2, sum, na.rm = TRUE)

      # From HF2 to Other states
      # HF2 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF2", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF2", "HFdeath", x] > 1] <- 1
      }

      # HF2 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF2", "Noncvdeath", x] > 1] <- 1
      }

      # HF2 to HF3
      a_transition_matrix["HF2", "HF3", ] <- 1 - apply(a_transition_matrix["HF2", , ], 2, sum, na.rm = TRUE)

      # From HF3 to Other states
      # HF3 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF3", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF3", "HFdeath", x] > 1] <- 1
      }

      # HF3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF3", "Noncvdeath", x] > 1] <- 1
      }

      # HF3 to HF4
      a_transition_matrix["HF3", "HF4", ] <- 1 - apply(a_transition_matrix["HF3", , ], 2, sum, na.rm = TRUE)

      # From HF4 to Other states
      # HF4 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF4", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF4", "HFdeath", x] > 1] <- 1
      }

      # HF4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF4", "Noncvdeath", x] > 1] <- 1
      }

      # HF4 to HF5
      a_transition_matrix["HF4", "HF5", ] <- 1 - apply(a_transition_matrix["HF4", , ], 2, sum, na.rm = TRUE)

      # From HF5 to Other states
      # HF5 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF5", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF5", "HFdeath", x] > 1] <- 1
      }

      # HF5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF5", "Noncvdeath", x] > 1] <- 1
      }

      # HF5 to HF6
      a_transition_matrix["HF5", "HF6", ] <- 1 - apply(a_transition_matrix["HF5", , ], 2, sum, na.rm = TRUE)

      # From HF6 to Other states
      # HF6 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF6", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF6", "HFdeath", x] > 1] <- 1
      }

      # HF6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF6", "Noncvdeath", x] > 1] <- 1
      }

      # HF6 to HF7
      a_transition_matrix["HF6", "HF7", ] <- 1 - apply(a_transition_matrix["HF6", , ], 2, sum, na.rm = TRUE)

      # From HF7 to Other states
      # HF7 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF7", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF7", "HFdeath", x] > 1] <- 1
      }

      # HF7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF7", "Noncvdeath", x] > 1] <- 1
      }

      # HF7 to HF8
      a_transition_matrix["HF7", "HF8", ] <- 1 - apply(a_transition_matrix["HF7", , ], 2, sum, na.rm = TRUE)

      # From HF8 to Other states
      # HF8 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF8", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF8", "HFdeath", x] > 1] <- 1
      }

      # HF8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF8", "Noncvdeath", x] > 1] <- 1
      }

      # HF8 to HF9
      a_transition_matrix["HF8", "HF9", ] <- 1 - apply(a_transition_matrix["HF8", , ], 2, sum, na.rm = TRUE)

      # From HF9 to Other states
      # HF9 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF9", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF9", "HFdeath", x] > 1] <- 1
      }

      # HF9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF9", "Noncvdeath", x] > 1] <- 1
      }

      # HF9 to HF10
      a_transition_matrix["HF9", "HF10", ] <- 1 - apply(a_transition_matrix["HF9", , ], 2, sum, na.rm = TRUE)

      # From HF10 to Other states
      # HF10 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF10", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF10", "HFdeath", x] > 1] <- 1
      }

      # HF10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF10", "Noncvdeath", x] > 1] <- 1
      }

      # HF10 to HF11
      a_transition_matrix["HF10", "HF11", ] <- 1 - apply(a_transition_matrix["HF10", , ], 2, sum, na.rm = TRUE)

      # From HF11 to Other states
      # HF11 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF11", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF11", "HFdeath", x] > 1] <- 1
      }

      # HF11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF11", "Noncvdeath", x] > 1] <- 1
      }

      # HF11 to HF12
      a_transition_matrix["HF11", "HF12", ] <- 1 - apply(a_transition_matrix["HF11", , ], 2, sum, na.rm = TRUE)

      # From HF12 to Other states
      # HF12 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF12", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["HF12", "HFdeath", x] > 1] <- 1
      }

      # HF12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["HF12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["HF12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix["HF12", "MI", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_MI")]
      a_transition_matrix["HF12", "HF", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_HF")]
      a_transition_matrix["HF12", "Stroke", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Stroke")]
      a_transition_matrix["HF12", "TIA", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_TIA")]
      a_transition_matrix["HF12", "Othercvd", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HF_Othercvd")]
      a_transition_matrix["HF12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix["HF12", , ], 2, sum, na.rm = TRUE)

      # From Stroke to Other states
      # Stroke to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke to Stroke2
      a_transition_matrix["Stroke", "Stroke2", ] <- 1 - apply(a_transition_matrix["Stroke", , ], 2, sum, na.rm = TRUE)

      # From Stroke2 to Other states
      # Stroke2 to Non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke2", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke2 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke2", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke2 to stroke3
      a_transition_matrix["Stroke2", "Stroke3", ] <- 1 - apply(a_transition_matrix["Stroke2", , ], 2, sum, na.rm = TRUE)

      # From Stroke3 to Other states
      # Stroke3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke3", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke3 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke3", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke3 to stroke4
      a_transition_matrix["Stroke3", "Stroke4", ] <- 1 - apply(a_transition_matrix["Stroke3", , ], 2, sum, na.rm = TRUE)

      # From Stroke4 to Other states
      # Stroke4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke4", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke4 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke4", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke4 to Stroke5
      a_transition_matrix["Stroke4", "Stroke5", ] <- 1 - apply(a_transition_matrix["Stroke4", , ], 2, sum, na.rm = TRUE)

      # From Stroke5 to Other states
      # Stroke5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke5", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke5 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke5", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke5 to Stroke6
      a_transition_matrix["Stroke5", "Stroke6", ] <- 1 - apply(a_transition_matrix["Stroke5", , ], 2, sum, na.rm = TRUE)

      # From Stroke6 to Other states
      # Stroke6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke6", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke6 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke6", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke6 to Stroke7
      a_transition_matrix["Stroke6", "Stroke7", ] <- 1 - apply(a_transition_matrix["Stroke6", , ], 2, sum, na.rm = TRUE)

      # From Stroke7 to Other states

      # Stroke7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke7", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke7 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke7", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke7 to Stroke8
      a_transition_matrix["Stroke7", "Stroke8", ] <- 1 - apply(a_transition_matrix["Stroke7", , ], 2, sum, na.rm = TRUE)

      # From Stroke8 to Other states
      # Stroke8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke8", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke8 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke8", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke8 to Stroke9
      a_transition_matrix["Stroke8", "Stroke9", ] <- 1 - apply(a_transition_matrix["Stroke8", , ], 2, sum, na.rm = TRUE)

      # From Stroke9 to Other states
      # Stroke9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke9", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke9 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke9", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke9 to Stroke10
      a_transition_matrix["Stroke9", "Stroke10", ] <- 1 - apply(a_transition_matrix["Stroke9", , ], 2, sum, na.rm = TRUE)

      # From Stroke10 to Other states

      # Stroke10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke10", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke10 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke10", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke10 to Stroke11
      a_transition_matrix["Stroke10", "Stroke11", ] <- 1 - apply(a_transition_matrix["Stroke10", , ], 2, sum, na.rm = TRUE)

      # From Stroke11 to Other states
      # Stroke11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke11", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke11 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke11", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke11 to Stroke12
      a_transition_matrix["Stroke11", "Stroke12", ] <- 1 - apply(a_transition_matrix["Stroke11", , ], 2, sum, na.rm = TRUE)

      # From Stroke12 to Other states
      # Stroke12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Stroke12", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke12 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["Stroke12", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      a_transition_matrix["Stroke12", "MI", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_MI")]
      a_transition_matrix["Stroke12", "HF", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_HF")]
      a_transition_matrix["Stroke12", "Stroke", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Stroke")]
      a_transition_matrix["Stroke12", "TIA", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_TIA")]
      a_transition_matrix["Stroke12", "Othercvd", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Stroke_Othercvd")]
      a_transition_matrix["Stroke12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix["Stroke12", , ], 2, sum, na.rm = TRUE)

      # From TIA to Other states
      # TIA to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA", "Noncvdeath", x] > 1] <- 1
      }

      # TIA to TIA2
      a_transition_matrix["TIA", "TIA2", ] <- 1 - apply(a_transition_matrix["TIA", , ], 2, sum, na.rm = TRUE)

      # From TIA2 to Other states
      # TIA2 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA2", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA2 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA2", "Noncvdeath", x] > 1] <- 1
      }

      # TIA2 to TIA3
      a_transition_matrix["TIA2", "TIA3", ] <- 1 - apply(a_transition_matrix["TIA2", , ], 2, sum, na.rm = TRUE)

      # From TIA3 to Other states
      # TIA3 to Stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA3", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA3", "Noncvdeath", x] > 1] <- 1
      }

      # TIA3 to TIa4
      a_transition_matrix["TIA3", "TIA4", ] <- 1 - apply(a_transition_matrix["TIA3", , ], 2, sum, na.rm = TRUE)

      # From TIA4 to Other states
      # TIA4 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA4", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA4", "Noncvdeath", x] > 1] <- 1
      }

      # TIA4 to TIA5
      a_transition_matrix["TIA4", "TIA5", ] <- 1 - apply(a_transition_matrix["TIA4", , ], 2, sum, na.rm = TRUE)

      # From TIA5 to Other states
      # TIA5 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA5", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA5", "Noncvdeath", x] > 1] <- 1
      }

      # TIA5 to TIA6
      a_transition_matrix["TIA5", "TIA6", ] <- 1 - apply(a_transition_matrix["TIA5", , ], 2, sum, na.rm = TRUE)

      # From TIA6 to Other states
      # TIA6 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA6", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA6", "Noncvdeath", x] > 1] <- 1
      }

      # TIA6 to TIA7
      a_transition_matrix["TIA6", "TIA7", ] <- 1 - apply(a_transition_matrix["TIA6", , ], 2, sum, na.rm = TRUE)

      # From TIA7 to Other states
      # TIA7 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA7", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA7  to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA7", "Noncvdeath", x] > 1] <- 1
      }

      # TIA7 to TIA8
      a_transition_matrix["TIA7", "TIA8", ] <- 1 - apply(a_transition_matrix["TIA7", , ], 2, sum, na.rm = TRUE)

      # From TIA8 to Other states
      # TIA8 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA8", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA8", "Noncvdeath", x] > 1] <- 1
      }

      # TIA8 to TIA9
      a_transition_matrix["TIA8", "TIA9", ] <- 1 - apply(a_transition_matrix["TIA8", , ], 2, sum, na.rm = TRUE)

      # From TIA9 to Other states
      # TIA9 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA9", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA9", "Noncvdeath", x] > 1] <- 1
      }

      # TIA9 to TIA10
      a_transition_matrix["TIA9", "TIA10", ] <- 1 - apply(a_transition_matrix["TIA9", , ], 2, sum, na.rm = TRUE)

      # From TIA10 to Other states
      # TIA10 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA10", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA10", "Noncvdeath", x] > 1] <- 1
      }

      # TIA10 to TIA11
      a_transition_matrix["TIA10", "TIA11", ] <- 1 - apply(a_transition_matrix["TIA10", , ], 2, sum, na.rm = TRUE)

      # From TIA11 to Other states
      # TIA11 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA11", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA11", "Noncvdeath", x] > 1] <- 1
      }

      # TIA11 to TIA12
      a_transition_matrix["TIA11", "TIA12", ] <- 1 - apply(a_transition_matrix["TIA11", , ], 2, sum, na.rm = TRUE)

      # From TIA12 to Other states
      # TIA12 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA12", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["TIA12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["TIA12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix["TIA12", "MI", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_MI")]
      a_transition_matrix["TIA12", "HF", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_HF")]
      a_transition_matrix["TIA12", "Stroke", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Stroke")]
      a_transition_matrix["TIA12", "TIA", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_TIA")]
      a_transition_matrix["TIA12", "Othercvd", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_TIA_Othercvd")]
      a_transition_matrix["TIA12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix["TIA12", , ], 2, sum, na.rm = TRUE)

      # From Othercvd to other states
      # Othercvd to non-cv death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd", "Othercvddeath", x] > 1] <- 1
      }
      # Othercvd to Othercvd2
      a_transition_matrix["Othercvd", "Othercvd2", ] <- 1 - apply(a_transition_matrix["Othercvd", , ], 2, sum, na.rm = TRUE)

      # From Othercvd2 to Other states
      # Othercvd2 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd2", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd2", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd2 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd2", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd2 to Othercvd3
      a_transition_matrix["Othercvd2", "Othercvd3", ] <- 1 - apply(a_transition_matrix["Othercvd2", , ], 2, sum, na.rm = TRUE)

      # From Othercvd3 to Other states
      # Othercvd3 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd3", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd3", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd3", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd3 to Othercvd4
      a_transition_matrix["Othercvd3", "Othercvd4", ] <- 1 - apply(a_transition_matrix["Othercvd3", , ], 2, sum, na.rm = TRUE)

      # From Othercvd4 to Other states
      # Othercvd4 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd4", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd4", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd4", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd4 to Othercvd5
      a_transition_matrix["Othercvd4", "Othercvd5", ] <- 1 - apply(a_transition_matrix["Othercvd4", , ], 2, sum, na.rm = TRUE)

      # From Othercvd5 to Other states
      # Othercvd5 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd5", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd5", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd5", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd5 to Othercvd6
      a_transition_matrix["Othercvd5", "Othercvd6", ] <- 1 - apply(a_transition_matrix["Othercvd5", , ], 2, sum, na.rm = TRUE)

      # From Othercvd6 to Other states
      # Othercvd6 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd6", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd6", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd6", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd6 to Othercvd7
      a_transition_matrix["Othercvd6", "Othercvd7", ] <- 1 - apply(a_transition_matrix["Othercvd6", , ], 2, sum, na.rm = TRUE)

      # From Othercvd7 to Other states
      # Othercvd7 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd7", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd7", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd7", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd7 to Othercvd8
      a_transition_matrix["Othercvd7", "Othercvd8", ] <- 1 - apply(a_transition_matrix["Othercvd7", , ], 2, sum, na.rm = TRUE)

      # From Othercvd8 to Other states
      # Othercvd8 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd8", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd8", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd8", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd8 to Othercvd9
      a_transition_matrix["Othercvd8", "Othercvd9", ] <- 1 - apply(a_transition_matrix["Othercvd8", , ], 2, sum, na.rm = TRUE)

      # From Othercvd9 to Other states
      # Othercvd9 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd9", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd9", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd9", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd9 to Othercvd10
      a_transition_matrix["Othercvd9", "Othercvd10", ] <- 1 - apply(a_transition_matrix["Othercvd9", , ], 2, sum, na.rm = TRUE)

      # From Othercvd10 to Other states
      # Othercvd10 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd10", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["TIA10", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd10", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd10 to Othercv11
      a_transition_matrix["Othercvd10", "Othercvd11", ] <- 1 - apply(a_transition_matrix["Othercvd10", , ], 2, sum, na.rm = TRUE)

      # From Othercvd11 to Other states
      # Othercvd11 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd11", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd11", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd11", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd11 to Othercvd12
      a_transition_matrix["Othercvd11", "Othercvd12", ] <- 1 - apply(a_transition_matrix["Othercvd11", , ], 2, sum, na.rm = TRUE)

      # From Othercvd12 to Other states
      # Othercvd12 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd12", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Othercvd12", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Othercvd12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Othercvd12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix["Othercvd12", "MI", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_MI")]
      a_transition_matrix["Othercvd12", "HF", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_HF")]
      a_transition_matrix["Othercvd12", "Stroke", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Stroke")]
      a_transition_matrix["Othercvd12", "TIA", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_TIA")]
      a_transition_matrix["Othercvd12", "Othercvd", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Othercvd_Othercvd")]
      a_transition_matrix["Othercvd12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix["Othercvd12", , ], 2, sum, na.rm = TRUE)


      # From Renal disease to Other states
      # Renal to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix["Renal", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"] * hr_renal
        a_transition_matrix[a_transition_matrix["Renal", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix["Renal", "Advrenal", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Renal_Advrenal")]
      a_transition_matrix["Renal", "Renal", ] <- 1 - apply(a_transition_matrix["Renal", , ], 2, sum, na.rm = TRUE)

      # From Chroniccvd to Other states
      # Chronic cvd to non-cv death
      for (x in 1:n_cycles) {
        a_transition_matrix["Chroniccvd", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"]
        a_transition_matrix[a_transition_matrix["Chroniccvd", "Noncvdeath", x] > 1] <- 1
      }

      # Chronic cvd to Chroniccvd death
      for (x in 1:n_cycles) {
        a_transition_matrix["Chroniccvd", "Chroniccvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_chroniccvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix[a_transition_matrix["Chroniccvd", "Chroniccvddeath", x] > 1] <- 1
      }

      a_transition_matrix["Chroniccvd", "MI", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_MI")]
      a_transition_matrix["Chroniccvd", "HF", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_HF")]
      a_transition_matrix["Chroniccvd", "Stroke", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_Stroke")]
      a_transition_matrix["Chroniccvd", "TIA", ] <- transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Chroniccvd_TIA")]
      a_transition_matrix["Chroniccvd", "Chroniccvd", ] <- 1 - apply(a_transition_matrix["Chroniccvd", , ], 2, sum, na.rm = TRUE)

      # From Advanced renal to Other states

      # Advanced renal to non-cv death
      for (x in 1:n_cycles) {
        a_transition_matrix["Advrenal", "Noncvdeath", x] <- (background_mortality$mnth_prob[x] * m_lhs_param_values[cb, "back_multiplier"] * hr_advrenal)
        a_transition_matrix[a_transition_matrix["Advrenal", "Noncvdeath", x] > 1] <- 1
      }
      a_transition_matrix["Advrenal", "Advrenal", ] <- 1 - apply(a_transition_matrix["Advrenal", , ], 2, sum, na.rm = TRUE)

      # From MI death to Death
      a_transition_matrix["MIdeath", "Death", ] <- 1

      # From HF death to Death
      a_transition_matrix["HFdeath", "Death", ] <- 1

      # From stroke death to Death
      a_transition_matrix["Strokedeath", "Death", ] <- 1

      # From Othercvd death to Death
      a_transition_matrix["Othercvddeath", "Death", ] <- 1

      # From non-CV death to Death
      a_transition_matrix["Noncvdeath", "Death", ] <- 1

      # From chroniccvd death to Death
      a_transition_matrix["Chroniccvddeath", "Death", ] <- 1

      # From Death to Other states
      a_transition_matrix["Death", "Death", ] <- 1

      # Replace NAs with 0 in the transition_matrix
      a_transition_matrix[is.na(a_transition_matrix)] <- 0

      # Standardizing the probabilities so that they sum up to 1
      for (x in 1:n_cycles) {
        a_transition_matrix[, , x] <- a_transition_matrix[, , x] / (apply(a_transition_matrix[, , x], MARGIN = 1, FUN = sum, na.rm = TRUE))
      }
      return(a_transition_matrix)
    })
  }

  # Calling the function for transition probabilities####
  m_transitionprobs <- prob_func(l_hrparams)


  # Multiplying the initial cohort with transition matrix####
  for (n in 1:n.i) {
    v.age.monthsstart <- v.age.monthsfrom60[n]
    a_state_membership[v.age.monthsstart, ] <- c(
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    for (i in (v.age.monthsstart + 1):n_cycles) {
      a_state_membership[i, ] <- a_state_membership[i - 1, ] %*% m_transitionprobs[, , i - 1]
      index.state <- sample(cond1.index, prob = a_state_membership[i, ], size = 1)
      m.prob[n, i - v.age.monthsstart + 1] <- cond1[index.state]
      a_state_membership[i, ] <- c(
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      )
      a_state_membership[i, index.state] <- 1
    }
  }

  # Replacing all NAs with Death
  m.prob1 <- m.prob
  m.prob1[is.na(m.prob1)] <- "Death"

  # Converting to dataframe
  m.prob2 <- as.data.frame(m.prob1)

  # Calculating how many patients died in cycles 1 - 169
  m_gof[cb, "n_Death_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:169)] == "Death", 1, any), ])

  # Calculating how many patients had CV death in cycles 1 - 168
  m_gof[cb, "n_Cvdeath_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:168)] == "MIdeath", 1, any), ]) +
    nrow(m.prob2[apply(m.prob2[c(1:168)] == "HFdeath", 1, any), ]) +
    nrow(m.prob2[apply(m.prob2[c(1:168)] == "Strokedeath", 1, any), ]) +
    nrow(m.prob2[apply(m.prob2[c(1:168)] == "Othercvddeath", 1, any), ]) +
    nrow(m.prob2[apply(m.prob2[c(1:168)] == "Chroniccvddeath", 1, any), ])

  # Calculating how many patients had non-CV death in cycles 1 - 168
  m_gof[cb, "n_Noncvdeath_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:168)] == "Noncvdeath", 1, any), ])

  # Calculating how many patients had stroke death in cycles 1 - 264
  m_gof[cb, "n_Strokedeath_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:264)] == "Strokedeath", 1, any), ])

  # Calculating how many patients had stroke in cycles 1- 60
  m_gof[cb, "n_Stroke_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:60)] == "Stroke", 1, any), ])

  # Calculating how many patients had TIA in cycles 1 - 60
  m_gof[cb, "n_TIA_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:60)] == "TIA", 1, any), ])

  # Calculating how many patients had MI in cycles 1 - 60
  m_gof[cb, "n_MI_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:60)] == "MI", 1, any), ])

  # Calculating how many patients had HF in cycles 1 - 60
  m_gof[cb, "n_HF_fit"] <- nrow(m.prob2[apply(m.prob2[c(1:60)] == "HF", 1, any), ])

  m_gof[cb, ]
}

stopCluster(numclusters)

# Calculating Poisson deviance for all calibration targets
df_lhscalibration_results1 <- as.data.frame(df_lhscalibration_results)

poisson_deviance <- function(obs, sim) {
  return(2 * ((obs * (log(obs / sim))) - (obs - sim)))
}

df_lhscalibration_results1$deviance_n_death <- mapply(poisson_deviance, df_lhscalibration_results1$n_Death, df_lhscalibration_results1$n_Death_fit)
df_lhscalibration_results1$deviance_n_cvdeath <- mapply(poisson_deviance, df_lhscalibration_results1$n_Cvdeath, df_lhscalibration_results1$n_Cvdeath_fit)
df_lhscalibration_results1$deviance_n_noncvdeath <- mapply(poisson_deviance, df_lhscalibration_results1$n_Noncvdeath, df_lhscalibration_results1$n_Noncvdeath_fit)
df_lhscalibration_results1$deviance_n_strokedeath <- mapply(poisson_deviance, df_lhscalibration_results1$n_Strokedeath, df_lhscalibration_results1$n_Strokedeath_fit)
df_lhscalibration_results1$deviance_n_mi <- mapply(poisson_deviance, df_lhscalibration_results1$n_MI, df_lhscalibration_results1$n_MI_fit)
df_lhscalibration_results1$deviance_n_stroke <- mapply(poisson_deviance, df_lhscalibration_results1$n_Stroke, df_lhscalibration_results1$n_Stroke_fit)
df_lhscalibration_results1$deviance_n_tia <- mapply(poisson_deviance, df_lhscalibration_results1$n_TIA, df_lhscalibration_results1$n_TIA_fit)
df_lhscalibration_results1$deviance_n_hf <- mapply(poisson_deviance, df_lhscalibration_results1$n_HF, df_lhscalibration_results1$n_HF_fit)

# Adding up all the deviance scores
df_lhscalibration_results1$deviance_all <- with(df_lhscalibration_results1, deviance_n_death + deviance_n_cvdeath + deviance_n_noncvdeath +
  deviance_n_strokedeath + deviance_n_mi + deviance_n_stroke + deviance_n_tia + deviance_n_hf)

# Reordering based on ascending order of total deviance score
df_lhscalibration_results1 <- df_lhscalibration_results1[order(df_lhscalibration_results1$deviance_all), ]

# Saving the df_lhscalibration_results1 R data
save(df_lhscalibration_results1, file = here("Outputs", paste0("df_lhscalibration_results1_", Sys.Date(), ".RData")))

# Saving the m_lhs_param_values R data
save(m_lhs_param_values, file = here("Outputs", paste0("m_lhs_param_values_", Sys.Date(), ".RData")))
