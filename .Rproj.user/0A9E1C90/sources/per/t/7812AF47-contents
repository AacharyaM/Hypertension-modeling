#######################################################################
# Filename: 02_Microsimulation_calibrate_nm.R
# Created by: Mahip
# Date created: 05/12/2025
# Purpose: Creating a Markov microsimulation model for
# patients >=60 years to calibrate the parameters
# 2nd round of calibration (Recalibration): Nelder-Mead algorithm
# Externally validating the model
#######################################################################

# Loading all the required libraries (installing them if they are not installed already)
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

# Load the individual libraries
loadpackages("here")
loadpackages("parallel")
loadpackages("dplyr")
loadpackages("tidyr")
loadpackages("ggplot2")
loadpackages("viridis")
loadpackages("ggrepel")

# Also installing readxl, truncnorm, readr, and dfoptim package but not calling the entire library and only using where needed
installonly <- function(packagename) {
  if (!requireNamespace(packagename, quietly = TRUE)) {
    install.packages(packagename)
  }
}

# Calling the function
installonly("readxl")
installonly("truncnorm")
installonly("readr")
installonly("dfoptim")

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


# Creating a vector of numbers corresponding to disease states
cond1.index <- c(1:71)

# Number of individuals
n.i <- 100000


# Creating a vector of n.i individuals with age drawn from a truncated
# distribution
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


####################################################################

# Recalibrating using Nelder-Mead algorithm######################

####################################################################

# The model was not calibrated for stroke deaths and HF incidence. Recalibrating the model using Nelder-Mead simplex method changing the transition
# probabilities to HF and mortality rate from stroke

# Loading the files from 1st round of calibration df_lhscalibration_results1 R data (dataframe where each row is a result of lhs calibration)
# Loading based on the latest date on the previously saved R objects
func_loadfiles <- function(Robject) {
  list_calib_filenames <- list.files(here("Outputs"), pattern = Robject)
  extract_dates <- gsub("[^0-9]", "", list_calib_filenames)
  latest_file <- list_calib_filenames[which.max(extract_dates)]
  load(here("Outputs", latest_file), envir = .GlobalEnv)
}

# Loading the files from 1st round of calibration df_lhscalibration_results1 R data (dataframe where each row is a result of lhs calibration)
func_loadfiles("df_lhscalibration_results")

# Loading the m_lhs_param_values R data (dataframe where each row corresponds to input values for the lhs calibration)
func_loadfiles("m_lhs_param_values")

# Identifying the iteration number that led to the best-fitting model
forindex.m_lhs_param_values <- rownames(df_lhscalibration_results1[1, ])
index.m_lhs_param_values <- readr::parse_number(gsub("result.", "", forindex.m_lhs_param_values))

# Extracting the input values for the best-fitting model
m_lhs_param_values1 <- m_lhs_param_values[index.m_lhs_param_values, ]

# Range of input parameters for changing in calibration
lowbound_nm <- c(0.7 * m_lhs_param_values1["p_HTN_HF"],
  0.7 * m_lhs_param_values1["p_MI_HF"],
  0.7 * m_lhs_param_values1["p_HF_HF"],
  0.7 * m_lhs_param_values1["p_Stroke_HF"],
  0.7 * m_lhs_param_values1["p_TIA_HF"],
  0.7 * m_lhs_param_values1["p_Othercvd_HF"],
  0.7 * m_lhs_param_values1["p_Chroniccvd_HF"],
  hr_stroke_multiplier = 1.2
)

upperbound_nm <- c(0.9 * m_lhs_param_values1["p_HTN_HF"],
  0.9 * m_lhs_param_values1["p_MI_HF"],
  0.9 * m_lhs_param_values1["p_HF_HF"],
  0.9 * m_lhs_param_values1["p_Stroke_HF"],
  0.9 * m_lhs_param_values1["p_TIA_HF"],
  0.9 * m_lhs_param_values1["p_Othercvd_HF"],
  0.9 * m_lhs_param_values1["p_Chroniccvd_HF"],
  hr_stroke_multiplier = 1.6
)

# Sample multiple random start values for the parameters to be calibrated using the
# Nelder-Mead algorithm
n_init_nm <- 10
n_param_nm <- 8
v_params_init_nm <- matrix(nrow = n_init_nm, ncol = n_param_nm)
for (i in 1:n_param_nm) {
  v_params_init_nm[, i] <- runif(n_init_nm, min = lowbound_nm[i], max = upperbound_nm[i])
}

v_hr_stroke_nm <- matrix(v_params_init_nm[, 8] * l_hrparams$hr_stroke * m_lhs_param_values1["hr_multiplier"], nrow = 10, ncol = 1)

# Combining the two
v_params_init1_nm <- cbind(v_params_init_nm, v_hr_stroke_nm)
colnames(v_params_init1_nm) <- c("p_HTN_HF", "p_MI_HF", "p_HF_HF", "p_Stroke_HF", "p_TIA_HF", "p_Othercvd_HF", "p_Chroniccvd_HF", "hr_stroke_multiplier", "hr_stroke")

# Removing the hr_stroke_multiplier
v_params_init2_nm <- v_params_init1_nm[, colnames(v_params_init1_nm) != "hr_stroke_multiplier"]
n_params_init2_nm <- ncol(v_params_init2_nm)
v_param_names_nm <- colnames(v_params_init2_nm)

# Creating a vector of the parameters that don't need calibrating
m_lhs_param_values2 <- m_lhs_param_values1[c(
  "back_multiplier", "p_HTN_MI", "p_HTN_Stroke", "p_HTN_TIA", "p_HTN_Othercvd", "p_MI_MI", "p_MI_Stroke", "p_MI_TIA", "p_MI_Othercvd", "p_HF_MI", "p_HF_Stroke",
  "p_HF_TIA", "p_HF_Othercvd", "p_Stroke_MI", "p_Stroke_Stroke", "p_Stroke_TIA", "p_Stroke_Othercvd", "p_TIA_MI", "p_TIA_Stroke", "p_TIA_TIA",
  "p_TIA_Othercvd", "p_Othercvd_MI", "p_Othercvd_Stroke", "p_Othercvd_TIA", "p_Othercvd_Othercvd", "p_Chroniccvd_MI", "p_Chroniccvd_Stroke", "p_Chroniccvd_TIA"
)]
# Converting to matrix
m_lhs_param_values3 <- matrix(m_lhs_param_values2, nrow = 1)

# Adding the parameters that were not used in the LHS calibration
m_noncvd_params <- c(p_HTN_Renal = transition_rates$mnth_prob[which(transition_rates$prob_type == "p_HTN_Renal")], p_Renal_Advrenal = transition_rates$mnth_prob[which(transition_rates$prob_type == "p_Renal_Advrenal")])

m_noncvd_params1 <- matrix(m_noncvd_params, nrow = 1)

# Adding the hazard ratios
m_lhs_param_values2a <- c(
  hr_mi = l_hrparams$hr_mi * m_lhs_param_values1["hr_multiplier"], hr_hf = l_hrparams$hr_hf * m_lhs_param_values1["hr_multiplier"], hr_tia = l_hrparams$hr_tia * m_lhs_param_values1["hr_multiplier"],
  hr_othercvd = l_hrparams$hr_othercvd * m_lhs_param_values1["hr_multiplier"], hr_chroniccvd = l_hrparams$hr_chroniccvd * m_lhs_param_values1["hr_multiplier"],
  hr_renal = l_hrparams$hr_renal, hr_advrenal = l_hrparams$hr_advrenal
)

m_lhs_param_values2b <- matrix(m_lhs_param_values2a, nrow = 1)

# Combining the three vectors (as a matrix)
m_lhs_param_values4 <- cbind(m_lhs_param_values3, m_noncvd_params1, m_lhs_param_values2b)

colnames(m_lhs_param_values4) <- c(
  "back_multiplier", "p_HTN_MI", "p_HTN_Stroke", "p_HTN_TIA", "p_HTN_Othercvd", "p_MI_MI", "p_MI_Stroke", "p_MI_TIA", "p_MI_Othercvd", "p_HF_MI", "p_HF_Stroke",
  "p_HF_TIA", "p_HF_Othercvd", "p_Stroke_MI", "p_Stroke_Stroke", "p_Stroke_TIA", "p_Stroke_Othercvd", "p_TIA_MI", "p_TIA_Stroke", "p_TIA_TIA",
  "p_TIA_Othercvd", "p_Othercvd_MI", "p_Othercvd_Stroke", "p_Othercvd_TIA", "p_Othercvd_Othercvd", "p_Chroniccvd_MI", "p_Chroniccvd_Stroke", "p_Chroniccvd_TIA",
  "p_HTN_Renal", "p_Renal_Advrenal", "hr_mi", "hr_hf", "hr_tia", "hr_othercvd", "hr_chroniccvd", "hr_renal", "hr_advrenal"
)

# Calibration targets
n_lhs_targets <- 8
v_lhs_targets <- c("n_Death", "n_MI", "n_HF", "n_Stroke", "n_TIA", "n_Cvdeath", "n_Noncvdeath", "n_Strokedeath")

m_gof_nm <- matrix(nrow = 1, ncol = n_lhs_targets * 2)

colnames(m_gof_nm)[9:16] <- paste0(v_lhs_targets, "_fit")

colnames(m_gof_nm)[1:8] <- v_lhs_targets

# Values for calibration targets
m_gof_nm[, "n_Death"] <- 45000 # 169 cycles
m_gof_nm[, "n_Cvdeath"] <- 24300 # 168 cycles
m_gof_nm[, "n_Noncvdeath"] <- 20700 # 168 cycles
m_gof_nm[, "n_Strokedeath"] <- 5000 # 264 cycles
m_gof_nm[, "n_MI"] <- 5054 # 60 cycles
m_gof_nm[, "n_HF"] <- 5483 # 60 cycles
m_gof_nm[, "n_Stroke"] <- 8200 # 60 cycles
m_gof_nm[, "n_TIA"] <- 4229 # 60 cycles


# Creating an empty transition matrix for Nelder-Mead
a_transition_matrix_nm <- array(NA_real_, dim = c(n_states, n_states, n_cycles), dimnames = list(
  from = cond1,
  to = cond1,
  cycle = 1:n_cycles
))

# Creating an empty array for the simulation cohort
a_state_membership_nm <- array(NA_real_,
  dim = c(n_cycles, n_states),
  dimnames = list(
    cycle = 1:n_cycles,
    state = cond1
  )
)

# Creating a matrix for states
m.prob_nm <- matrix(
  nrow = n.i, ncol = n_cycles,
  dimnames = list(
    paste("ind", 1:n.i, sep = " "),
    paste("cycle", 1:n_cycles, sep = " ")
  )
)

# Initializing the first cycle as HTN state
m.prob_nm[, 1] <- "HTN"

# Creating common random numbers
crn <- sample.int(1000000, 100000)

# Function for model and goodness of fit
prob_func_nm <- function(params) {
  with(as.list(params), {
    # Extracting the individual probabilities used in NM algorithm
    p_HTN_HF <- params[1]
    p_MI_HF <- params[2]
    p_HF_HF <- params[3]
    p_Stroke_HF <- params[4]
    p_TIA_HF <- params[5]
    p_Othercvd_HF <- params[6]
    p_Chroniccvd_HF <- params[7]
    hr_stroke <- params[8]

    ## Filling the transition probability values##

    with(as.list(as.data.frame(m_lhs_param_values4)), {
      # From HTN to Death across all cycles
      # HTN to non-cardiovascular death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HTN", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HTN", "Noncvdeath", x] > 1] <- 1
      }

      # From HTN to Other states
      a_transition_matrix_nm["HTN", "MI", ] <- p_HTN_MI
      a_transition_matrix_nm["HTN", "HF", ] <- p_HTN_HF
      a_transition_matrix_nm["HTN", "Stroke", ] <- p_HTN_Stroke
      a_transition_matrix_nm["HTN", "TIA", ] <- p_HTN_TIA
      a_transition_matrix_nm["HTN", "Renal", ] <- p_HTN_Renal
      a_transition_matrix_nm["HTN", "Othercvd", ] <- p_HTN_Othercvd
      a_transition_matrix_nm["HTN", "HTN", ] <- 1 - apply(a_transition_matrix_nm["HTN", , ], 2, sum, na.rm = TRUE)

      # From MI to all other states
      # MI to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI", "Noncvdeath", x] > 1] <- 1
      }

      # MI to MI2
      a_transition_matrix_nm["MI", "MI2", ] <- 1 - apply(a_transition_matrix_nm["MI", , ], 2, sum, na.rm = TRUE)

      # From MI2 to Other states
      # MI2 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI2", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI2 to Non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI2", "Noncvdeath", x] > 1] <- 1
      }

      # MI2 to MI3
      a_transition_matrix_nm["MI2", "MI3", ] <- 1 - apply(a_transition_matrix_nm["MI2", , ], 2, sum, na.rm = TRUE)

      # From MI3 to Other states
      # MI3 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI3", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI3", "Noncvdeath", x] > 1] <- 1
      }

      # MI3 to MI4
      a_transition_matrix_nm["MI3", "MI4", ] <- 1 - apply(a_transition_matrix_nm["MI3", , ], 2, sum, na.rm = TRUE)

      # From MI4 to Other states
      # MI4 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI4", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI4", "Noncvdeath", x] > 1] <- 1
      }

      # MI4 to MI5
      a_transition_matrix_nm["MI4", "MI5", ] <- 1 - apply(a_transition_matrix_nm["MI4", , ], 2, sum, na.rm = TRUE)

      # From MI5 to Other states
      # MI5 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI5", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI5", "Noncvdeath", x] > 1] <- 1
      }

      # MI5 to MI6
      a_transition_matrix_nm["MI5", "MI6", ] <- 1 - apply(a_transition_matrix_nm["MI5", , ], 2, sum, na.rm = TRUE)

      # From MI6 to Other states
      # MI6 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI6", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI6", "Noncvdeath", x] > 1] <- 1
      }

      # MI6 to MI7
      a_transition_matrix_nm["MI6", "MI7", ] <- 1 - apply(a_transition_matrix_nm["MI6", , ], 2, sum, na.rm = TRUE)

      # From MI7 to Other states
      # MI7 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI7", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI7", "Noncvdeath", x] > 1] <- 1
      }

      # MI7 to MI8
      a_transition_matrix_nm["MI7", "MI8", ] <- 1 - apply(a_transition_matrix_nm["MI7", , ], 2, sum, na.rm = TRUE)

      # From MI8 to Other states
      # MI8 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI8", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI8", "Noncvdeath", x] > 1] <- 1
      }

      # MI8 to MI9
      a_transition_matrix_nm["MI8", "MI9", ] <- 1 - apply(a_transition_matrix_nm["MI8", , ], 2, sum, na.rm = TRUE)

      # From MI9 to Other states
      # MI9 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI9", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI9", "Noncvdeath", x] > 1] <- 1
      }

      # MI9 to MI10
      a_transition_matrix_nm["MI9", "MI10", ] <- 1 - apply(a_transition_matrix_nm["MI9", , ], 2, sum, na.rm = TRUE)

      # From MI10 to Other states
      # MI10 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI10", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI10", "Noncvdeath", x] > 1] <- 1
      }

      # MI10 to MI11
      a_transition_matrix_nm["MI10", "MI11", ] <- 1 - apply(a_transition_matrix_nm["MI10", , ], 2, sum, na.rm = TRUE)

      # From MI11 to Other states
      # MI11 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI11", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # MI11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI11", "Noncvdeath", x] > 1] <- 1
      }

      # MI11 to MI12
      a_transition_matrix_nm["MI11", "MI12", ] <- 1 - apply(a_transition_matrix_nm["MI11", , ], 2, sum, na.rm = TRUE)

      # From MI12 to Other states
      # MI12 to MI death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI12", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["MI12", "MIdeath", x] > 1] <- 1
      }

      # MI12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["MI12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["MI12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix_nm["MI12", "MI", ] <- p_MI_MI
      a_transition_matrix_nm["MI12", "HF", ] <- p_MI_HF
      a_transition_matrix_nm["MI12", "Stroke", ] <- p_MI_Stroke
      a_transition_matrix_nm["MI12", "TIA", ] <- p_MI_TIA
      a_transition_matrix_nm["MI12", "Othercvd", ] <- p_MI_Othercvd
      a_transition_matrix_nm["MI12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_nm["MI12", , ], 2, sum, na.rm = TRUE)

      # From HF to Other states
      # HF to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF", "HFdeath", x] > 1] <- 1
      }

      # HF to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF", "Noncvdeath", x] > 1] <- 1
      }

      # HF to HF2
      a_transition_matrix_nm["HF", "HF2", ] <- 1 - apply(a_transition_matrix_nm["HF", , ], 2, sum, na.rm = TRUE)

      # From HF2 to Other states
      # HF2 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF2", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF2", "HFdeath", x] > 1] <- 1
      }

      # HF2 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF2", "Noncvdeath", x] > 1] <- 1
      }

      # HF2 to HF3
      a_transition_matrix_nm["HF2", "HF3", ] <- 1 - apply(a_transition_matrix_nm["HF2", , ], 2, sum, na.rm = TRUE)

      # From HF3 to Other states
      # HF3 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF3", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF3", "HFdeath", x] > 1] <- 1
      }

      # HF3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF3", "Noncvdeath", x] > 1] <- 1
      }

      # HF3 to HF4
      a_transition_matrix_nm["HF3", "HF4", ] <- 1 - apply(a_transition_matrix_nm["HF3", , ], 2, sum, na.rm = TRUE)

      # From HF4 to Other states
      # HF4 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF4", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF4", "HFdeath", x] > 1] <- 1
      }

      # HF4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF4", "Noncvdeath", x] > 1] <- 1
      }

      # HF4 to HF5
      a_transition_matrix_nm["HF4", "HF5", ] <- 1 - apply(a_transition_matrix_nm["HF4", , ], 2, sum, na.rm = TRUE)

      # From HF5 to Other states
      # HF5 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF5", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF5", "HFdeath", x] > 1] <- 1
      }

      # HF5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF5", "Noncvdeath", x] > 1] <- 1
      }

      # HF5 to HF6
      a_transition_matrix_nm["HF5", "HF6", ] <- 1 - apply(a_transition_matrix_nm["HF5", , ], 2, sum, na.rm = TRUE)

      # From HF6 to Other states
      # HF6 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF6", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF6", "HFdeath", x] > 1] <- 1
      }

      # HF6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF6", "Noncvdeath", x] > 1] <- 1
      }

      # HF6 to HF7
      a_transition_matrix_nm["HF6", "HF7", ] <- 1 - apply(a_transition_matrix_nm["HF6", , ], 2, sum, na.rm = TRUE)

      # From HF7 to Other states
      # HF7 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF7", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF7", "HFdeath", x] > 1] <- 1
      }

      # HF7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF7", "Noncvdeath", x] > 1] <- 1
      }

      # HF7 to HF8
      a_transition_matrix_nm["HF7", "HF8", ] <- 1 - apply(a_transition_matrix_nm["HF7", , ], 2, sum, na.rm = TRUE)

      # From HF8 to Other states
      # HF8 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF8", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF8", "HFdeath", x] > 1] <- 1
      }

      # HF8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF8", "Noncvdeath", x] > 1] <- 1
      }

      # HF8 to HF9
      a_transition_matrix_nm["HF8", "HF9", ] <- 1 - apply(a_transition_matrix_nm["HF8", , ], 2, sum, na.rm = TRUE)

      # From HF9 to Other states
      # HF9 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF9", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF9", "HFdeath", x] > 1] <- 1
      }

      # HF9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF9", "Noncvdeath", x] > 1] <- 1
      }

      # HF9 to HF10
      a_transition_matrix_nm["HF9", "HF10", ] <- 1 - apply(a_transition_matrix_nm["HF9", , ], 2, sum, na.rm = TRUE)

      # From HF10 to Other states
      # HF10 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF10", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF10", "HFdeath", x] > 1] <- 1
      }

      # HF10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF10", "Noncvdeath", x] > 1] <- 1
      }

      # HF10 to HF11
      a_transition_matrix_nm["HF10", "HF11", ] <- 1 - apply(a_transition_matrix_nm["HF10", , ], 2, sum, na.rm = TRUE)

      # From HF11 to Other states
      # HF11 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF11", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF11", "HFdeath", x] > 1] <- 1
      }

      # HF11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF11", "Noncvdeath", x] > 1] <- 1
      }

      # HF11 to HF12
      a_transition_matrix_nm["HF11", "HF12", ] <- 1 - apply(a_transition_matrix_nm["HF11", , ], 2, sum, na.rm = TRUE)

      # From HF12 to Other states
      # HF12 to HF death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF12", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["HF12", "HFdeath", x] > 1] <- 1
      }

      # HF12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["HF12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["HF12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix_nm["HF12", "MI", ] <- p_HF_MI
      a_transition_matrix_nm["HF12", "HF", ] <- p_HF_HF
      a_transition_matrix_nm["HF12", "Stroke", ] <- p_HF_Stroke
      a_transition_matrix_nm["HF12", "TIA", ] <- p_HF_TIA
      a_transition_matrix_nm["HF12", "Othercvd", ] <- p_HF_Othercvd
      a_transition_matrix_nm["HF12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_nm["HF12", , ], 2, sum, na.rm = TRUE)

      # From Stroke to Other states
      # Stroke to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke to Stroke2
      a_transition_matrix_nm["Stroke", "Stroke2", ] <- 1 - apply(a_transition_matrix_nm["Stroke", , ], 2, sum, na.rm = TRUE)

      # From Stroke2 to Other states
      # Stroke2 to Non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke2", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke2 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke2", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke2 to stroke3
      a_transition_matrix_nm["Stroke2", "Stroke3", ] <- 1 - apply(a_transition_matrix_nm["Stroke2", , ], 2, sum, na.rm = TRUE)

      # From Stroke3 to Other states
      # Stroke3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke3", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke3 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke3", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke3 to stroke4
      a_transition_matrix_nm["Stroke3", "Stroke4", ] <- 1 - apply(a_transition_matrix_nm["Stroke3", , ], 2, sum, na.rm = TRUE)

      # From Stroke4 to Other states
      # Stroke4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke4", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke4 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke4", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke4 to Stroke5
      a_transition_matrix_nm["Stroke4", "Stroke5", ] <- 1 - apply(a_transition_matrix_nm["Stroke4", , ], 2, sum, na.rm = TRUE)

      # From Stroke5 to Other states
      # Stroke5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke5", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke5 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke5", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke5 to Stroke6
      a_transition_matrix_nm["Stroke5", "Stroke6", ] <- 1 - apply(a_transition_matrix_nm["Stroke5", , ], 2, sum, na.rm = TRUE)

      # From Stroke6 to Other states
      # Stroke6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke6", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke6 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke6", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke6 to Stroke7
      a_transition_matrix_nm["Stroke6", "Stroke7", ] <- 1 - apply(a_transition_matrix_nm["Stroke6", , ], 2, sum, na.rm = TRUE)

      # From Stroke7 to Other states

      # Stroke7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke7", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke7 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke7", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke7 to Stroke8
      a_transition_matrix_nm["Stroke7", "Stroke8", ] <- 1 - apply(a_transition_matrix_nm["Stroke7", , ], 2, sum, na.rm = TRUE)

      # From Stroke8 to Other states
      # Stroke8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke8", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke8 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke8", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke8 to Stroke9
      a_transition_matrix_nm["Stroke8", "Stroke9", ] <- 1 - apply(a_transition_matrix_nm["Stroke8", , ], 2, sum, na.rm = TRUE)

      # From Stroke9 to Other states
      # Stroke9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke9", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke9 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke9", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke9 to Stroke10
      a_transition_matrix_nm["Stroke9", "Stroke10", ] <- 1 - apply(a_transition_matrix_nm["Stroke9", , ], 2, sum, na.rm = TRUE)

      # From Stroke10 to Other states

      # Stroke10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke10", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke10 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke10", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke10 to Stroke11
      a_transition_matrix_nm["Stroke10", "Stroke11", ] <- 1 - apply(a_transition_matrix_nm["Stroke10", , ], 2, sum, na.rm = TRUE)

      # From Stroke11 to Other states
      # Stroke11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke11", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke11 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke11", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # Stroke11 to Stroke12
      a_transition_matrix_nm["Stroke11", "Stroke12", ] <- 1 - apply(a_transition_matrix_nm["Stroke11", , ], 2, sum, na.rm = TRUE)

      # From Stroke12 to Other states
      # Stroke12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Stroke12", "Noncvdeath", x] > 1] <- 1
      }

      # Stroke12 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Stroke12", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      a_transition_matrix_nm["Stroke12", "MI", ] <- p_Stroke_MI
      a_transition_matrix_nm["Stroke12", "HF", ] <- p_Stroke_HF
      a_transition_matrix_nm["Stroke12", "Stroke", ] <- p_Stroke_Stroke
      a_transition_matrix_nm["Stroke12", "TIA", ] <- p_Stroke_TIA
      a_transition_matrix_nm["Stroke12", "Othercvd", ] <- p_Stroke_Othercvd
      a_transition_matrix_nm["Stroke12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_nm["Stroke12", , ], 2, sum, na.rm = TRUE)

      # From TIA to Other states
      # TIA to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA", "Noncvdeath", x] > 1] <- 1
      }

      # TIA to TIA2
      a_transition_matrix_nm["TIA", "TIA2", ] <- 1 - apply(a_transition_matrix_nm["TIA", , ], 2, sum, na.rm = TRUE)

      # From TIA2 to Other states
      # TIA2 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA2", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA2 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA2", "Noncvdeath", x] > 1] <- 1
      }

      # TIA2 to TIA3
      a_transition_matrix_nm["TIA2", "TIA3", ] <- 1 - apply(a_transition_matrix_nm["TIA2", , ], 2, sum, na.rm = TRUE)

      # From TIA3 to Other states
      # TIA3 to Stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA3", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA3", "Noncvdeath", x] > 1] <- 1
      }

      # TIA3 to TIa4
      a_transition_matrix_nm["TIA3", "TIA4", ] <- 1 - apply(a_transition_matrix_nm["TIA3", , ], 2, sum, na.rm = TRUE)

      # From TIA4 to Other states
      # TIA4 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA4", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA4", "Noncvdeath", x] > 1] <- 1
      }

      # TIA4 to TIA5
      a_transition_matrix_nm["TIA4", "TIA5", ] <- 1 - apply(a_transition_matrix_nm["TIA4", , ], 2, sum, na.rm = TRUE)

      # From TIA5 to Other states
      # TIA5 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA5", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA5", "Noncvdeath", x] > 1] <- 1
      }

      # TIA5 to TIA6
      a_transition_matrix_nm["TIA5", "TIA6", ] <- 1 - apply(a_transition_matrix_nm["TIA5", , ], 2, sum, na.rm = TRUE)

      # From TIA6 to Other states
      # TIA6 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA6", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA6", "Noncvdeath", x] > 1] <- 1
      }

      # TIA6 to TIA7
      a_transition_matrix_nm["TIA6", "TIA7", ] <- 1 - apply(a_transition_matrix_nm["TIA6", , ], 2, sum, na.rm = TRUE)

      # From TIA7 to Other states
      # TIA7 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA7", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA7  to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA7", "Noncvdeath", x] > 1] <- 1
      }

      # TIA7 to TIA8
      a_transition_matrix_nm["TIA7", "TIA8", ] <- 1 - apply(a_transition_matrix_nm["TIA7", , ], 2, sum, na.rm = TRUE)

      # From TIA8 to Other states
      # TIA8 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA8", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA8", "Noncvdeath", x] > 1] <- 1
      }

      # TIA8 to TIA9
      a_transition_matrix_nm["TIA8", "TIA9", ] <- 1 - apply(a_transition_matrix_nm["TIA8", , ], 2, sum, na.rm = TRUE)

      # From TIA9 to Other states
      # TIA9 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA9", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA9", "Noncvdeath", x] > 1] <- 1
      }

      # TIA9 to TIA10
      a_transition_matrix_nm["TIA9", "TIA10", ] <- 1 - apply(a_transition_matrix_nm["TIA9", , ], 2, sum, na.rm = TRUE)

      # From TIA10 to Other states
      # TIA10 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA10", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA10", "Noncvdeath", x] > 1] <- 1
      }

      # TIA10 to TIA11
      a_transition_matrix_nm["TIA10", "TIA11", ] <- 1 - apply(a_transition_matrix_nm["TIA10", , ], 2, sum, na.rm = TRUE)

      # From TIA11 to Other states
      # TIA11 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA11", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA11", "Noncvdeath", x] > 1] <- 1
      }

      # TIA11 to TIA12
      a_transition_matrix_nm["TIA11", "TIA12", ] <- 1 - apply(a_transition_matrix_nm["TIA11", , ], 2, sum, na.rm = TRUE)

      # From TIA12 to Other states
      # TIA12 to stroke death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA12", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      }

      # TIA12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["TIA12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["TIA12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix_nm["TIA12", "MI", ] <- p_TIA_MI
      a_transition_matrix_nm["TIA12", "HF", ] <- p_TIA_HF
      a_transition_matrix_nm["TIA12", "Stroke", ] <- p_TIA_Stroke
      a_transition_matrix_nm["TIA12", "TIA", ] <- p_TIA_TIA
      a_transition_matrix_nm["TIA12", "Othercvd", ] <- p_TIA_Othercvd
      a_transition_matrix_nm["TIA12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_nm["TIA12", , ], 2, sum, na.rm = TRUE)

      # From Othercvd to other states
      # Othercvd to non-cv death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd", "Othercvddeath", x] > 1] <- 1
      }
      # Othercvd to Othercvd2
      a_transition_matrix_nm["Othercvd", "Othercvd2", ] <- 1 - apply(a_transition_matrix_nm["Othercvd", , ], 2, sum, na.rm = TRUE)

      # From Othercvd2 to Other states
      # Othercvd2 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd2", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd2", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd2 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd2", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd2 to Othercvd3
      a_transition_matrix_nm["Othercvd2", "Othercvd3", ] <- 1 - apply(a_transition_matrix_nm["Othercvd2", , ], 2, sum, na.rm = TRUE)

      # From Othercvd3 to Other states
      # Othercvd3 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd3", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd3", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd3 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd3", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd3 to Othercvd4
      a_transition_matrix_nm["Othercvd3", "Othercvd4", ] <- 1 - apply(a_transition_matrix_nm["Othercvd3", , ], 2, sum, na.rm = TRUE)

      # From Othercvd4 to Other states
      # Othercvd4 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd4", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd4", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd4 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd4", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd4 to Othercvd5
      a_transition_matrix_nm["Othercvd4", "Othercvd5", ] <- 1 - apply(a_transition_matrix_nm["Othercvd4", , ], 2, sum, na.rm = TRUE)

      # From Othercvd5 to Other states
      # Othercvd5 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd5", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd5", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd5 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd5", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd5 to Othercvd6
      a_transition_matrix_nm["Othercvd5", "Othercvd6", ] <- 1 - apply(a_transition_matrix_nm["Othercvd5", , ], 2, sum, na.rm = TRUE)

      # From Othercvd6 to Other states
      # Othercvd6 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd6", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd6", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd6 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd6", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd6 to Othercvd7
      a_transition_matrix_nm["Othercvd6", "Othercvd7", ] <- 1 - apply(a_transition_matrix_nm["Othercvd6", , ], 2, sum, na.rm = TRUE)

      # From Othercvd7 to Other states
      # Othercvd7 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd7", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd7", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd7 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd7", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd7 to Othercvd8
      a_transition_matrix_nm["Othercvd7", "Othercvd8", ] <- 1 - apply(a_transition_matrix_nm["Othercvd7", , ], 2, sum, na.rm = TRUE)

      # From Othercvd8 to Other states
      # Othercvd8 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd8", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd8", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd8 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd8", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd8 to Othercvd9
      a_transition_matrix_nm["Othercvd8", "Othercvd9", ] <- 1 - apply(a_transition_matrix_nm["Othercvd8", , ], 2, sum, na.rm = TRUE)

      # From Othercvd9 to Other states
      # Othercvd9 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd9", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd9", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd9 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd9", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd9 to Othercvd10
      a_transition_matrix_nm["Othercvd9", "Othercvd10", ] <- 1 - apply(a_transition_matrix_nm["Othercvd9", , ], 2, sum, na.rm = TRUE)

      # From Othercvd10 to Other states
      # Othercvd10 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd10", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["TIA10", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd10 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd10", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd10 to Othercvd11
      a_transition_matrix_nm["Othercvd10", "Othercvd11", ] <- 1 - apply(a_transition_matrix_nm["Othercvd10", , ], 2, sum, na.rm = TRUE)

      # From Othercvd11 to Other states
      # Othercvd11 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd11", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd11", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd11 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd11", "Noncvdeath", x] > 1] <- 1
      }

      # Othercvd11 to Othercvd12
      a_transition_matrix_nm["Othercvd11", "Othercvd12", ] <- 1 - apply(a_transition_matrix_nm["Othercvd11", , ], 2, sum, na.rm = TRUE)

      # From Othercvd12 to Other states
      # Othercvd12 to Othercvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd12", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd12", "Othercvddeath", x] > 1] <- 1
      }

      # Othercvd12 to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Othercvd12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Othercvd12", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix_nm["Othercvd12", "MI", ] <- p_Othercvd_MI
      a_transition_matrix_nm["Othercvd12", "HF", ] <- p_Othercvd_HF
      a_transition_matrix_nm["Othercvd12", "Stroke", ] <- p_Othercvd_Stroke
      a_transition_matrix_nm["Othercvd12", "TIA", ] <- p_Othercvd_TIA
      a_transition_matrix_nm["Othercvd12", "Othercvd", ] <- p_Othercvd_Othercvd
      a_transition_matrix_nm["Othercvd12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_nm["Othercvd12", , ], 2, sum, na.rm = TRUE)


      # From Renal disease to Other states
      # Renal to non-CV death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Renal", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier * hr_renal
        a_transition_matrix_nm[a_transition_matrix_nm["Renal", "Noncvdeath", x] > 1] <- 1
      }

      a_transition_matrix_nm["Renal", "Advrenal", ] <- p_Renal_Advrenal
      a_transition_matrix_nm["Renal", "Renal", ] <- 1 - apply(a_transition_matrix_nm["Renal", , ], 2, sum, na.rm = TRUE)

      # From Chroniccvd to Other states
      # Chronic cvd to non-cv death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Chroniccvd", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
        a_transition_matrix_nm[a_transition_matrix_nm["Chroniccvd", "Noncvdeath", x] > 1] <- 1
      }

      # Chronic cvd to Chroniccvd death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Chroniccvd", "Chroniccvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_chroniccvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
        a_transition_matrix_nm[a_transition_matrix_nm["Chroniccvd", "Chroniccvddeath", x] > 1] <- 1
      }

      a_transition_matrix_nm["Chroniccvd", "MI", ] <- p_Chroniccvd_MI
      a_transition_matrix_nm["Chroniccvd", "HF", ] <- p_Chroniccvd_HF
      a_transition_matrix_nm["Chroniccvd", "Stroke", ] <- p_Chroniccvd_Stroke
      a_transition_matrix_nm["Chroniccvd", "TIA", ] <- p_Chroniccvd_TIA
      a_transition_matrix_nm["Chroniccvd", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_nm["Chroniccvd", , ], 2, sum, na.rm = TRUE)

      # From Advanced renal to Other states

      # Advanced renal to non-cv death
      for (x in 1:n_cycles) {
        a_transition_matrix_nm["Advrenal", "Noncvdeath", x] <- (background_mortality$mnth_prob[x] * back_multiplier * hr_advrenal)
        a_transition_matrix_nm[a_transition_matrix_nm["Advrenal", "Noncvdeath", x] > 1] <- 1
      }
      a_transition_matrix_nm["Advrenal", "Advrenal", ] <- 1 - apply(a_transition_matrix_nm["Advrenal", , ], 2, sum, na.rm = TRUE)

      # From MI death to Death
      a_transition_matrix_nm["MIdeath", "Death", ] <- 1

      # From HF death to Death
      a_transition_matrix_nm["HFdeath", "Death", ] <- 1

      # From stroke death to Death
      a_transition_matrix_nm["Strokedeath", "Death", ] <- 1

      # From Othercvd death to Death
      a_transition_matrix_nm["Othercvddeath", "Death", ] <- 1

      # From non-CV death to Death
      a_transition_matrix_nm["Noncvdeath", "Death", ] <- 1

      # From chroniccvd death to Death
      a_transition_matrix_nm["Chroniccvddeath", "Death", ] <- 1

      # From Death to Other states
      a_transition_matrix_nm["Death", "Death", ] <- 1

      # Replace NAs with 0 in the transition_matrix
      a_transition_matrix_nm[is.na(a_transition_matrix_nm)] <- 0

      # Standardizing the probabilities so that they sum up to 1
      for (x in 1:n_cycles) {
        a_transition_matrix_nm[, , x] <- a_transition_matrix_nm[, , x] / (apply(a_transition_matrix_nm[, , x], MARGIN = 1, FUN = sum, na.rm = TRUE))
      }
      return(a_transition_matrix_nm)
    })
  })
}

# Creating the function for microsimulation model with common random numbers####
microsim_nm <- function(params) {
  # Calling the function for transition probability matrix
  m_transitionprobs <- prob_func_nm(params)

  for (n in 1:n.i) {
    set.seed(crn[n])
    v.age.monthsstart <- v.age.monthsfrom60[n]
    a_state_membership_nm[v.age.monthsstart, ] <- c(
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    for (i in (v.age.monthsstart + 1):n_cycles) {
      a_state_membership_nm[i, ] <- a_state_membership_nm[i - 1, ] %*% m_transitionprobs[, , i - 1]
      index.state <- sample(cond1.index, prob = a_state_membership_nm[i, ], size = 1)
      m.prob_nm[n, i - v.age.monthsstart + 1] <- cond1[index.state]
      a_state_membership_nm[i, ] <- c(
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      )
      a_state_membership_nm[i, index.state] <- 1
    }
  }

  # Replacing all NAs with Death
  m.prob1_nm <- m.prob_nm
  m.prob1_nm[is.na(m.prob1_nm)] <- "Death"

  # Converting to dataframe
  m.prob2_nm <- as.data.frame(m.prob1_nm)

  # Calculating how many patients died in cycles 1 - 169
  m_gof_nm[, "n_Death_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:169)] == "Death", 1, any), ])

  # Calculating how many patients had CV death in cycles 1 - 168
  m_gof_nm[, "n_Cvdeath_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:168)] == "MIdeath", 1, any), ]) +
    nrow(m.prob2_nm[apply(m.prob2_nm[c(1:168)] == "HFdeath", 1, any), ]) +
    nrow(m.prob2_nm[apply(m.prob2_nm[c(1:168)] == "Strokedeath", 1, any), ]) +
    nrow(m.prob2_nm[apply(m.prob2_nm[c(1:168)] == "Othercvddeath", 1, any), ]) +
    nrow(m.prob2_nm[apply(m.prob2_nm[c(1:168)] == "Chroniccvddeath", 1, any), ])

  # Calculating how many patients had non-CV death in cycles 1 - 168
  m_gof_nm[, "n_Noncvdeath_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:168)] == "Noncvdeath", 1, any), ])

  # Calculating how many patients had stroke death in cycles 1 - 264
  m_gof_nm[, "n_Strokedeath_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:264)] == "Strokedeath", 1, any), ])

  # Calculating how many patients had stroke in cycles 1- 60
  m_gof_nm[, "n_Stroke_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:60)] == "Stroke", 1, any), ])

  # Calculating how many patients had TIA in cycles 1 - 60
  m_gof_nm[, "n_TIA_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:60)] == "TIA", 1, any), ])

  # Calculating how many patients had MI in cycles 1 - 60
  m_gof_nm[, "n_MI_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:60)] == "MI", 1, any), ])

  # Calculating how many patients had HF in cycles 1 - 60
  m_gof_nm[, "n_HF_fit"] <- nrow(m.prob2_nm[apply(m.prob2_nm[c(1:60)] == "HF", 1, any), ])

  # Calculating Poisson deviance for all calibration targets
  poisson_deviance <- function(obs, sim) {
    return(2 * ((obs * (log(obs / sim))) - (obs - sim)))
  }

  m_gof_nm1 <- as.data.frame(m_gof_nm)
  m_gof_nm1$deviance_n_death <- mapply(poisson_deviance, m_gof_nm1$n_Death, m_gof_nm1$n_Death_fit)
  m_gof_nm1$deviance_n_cvdeath <- mapply(poisson_deviance, m_gof_nm1$n_Cvdeath, m_gof_nm1$n_Cvdeath_fit)
  m_gof_nm1$deviance_n_noncvdeath <- mapply(poisson_deviance, m_gof_nm1$n_Noncvdeath, m_gof_nm1$n_Noncvdeath_fit)
  m_gof_nm1$deviance_n_mi <- mapply(poisson_deviance, m_gof_nm1$n_MI, m_gof_nm1$n_MI_fit)
  m_gof_nm1$deviance_n_stroke <- mapply(poisson_deviance, m_gof_nm1$n_Stroke, m_gof_nm1$n_Stroke_fit)
  m_gof_nm1$deviance_n_tia <- mapply(poisson_deviance, m_gof_nm1$n_TIA, m_gof_nm1$n_TIA_fit)
  m_gof_nm1$deviance_n_hf <- mapply(poisson_deviance, m_gof_nm1$n_HF, m_gof_nm1$n_HF_fit)
  m_gof_nm1$deviance_n_strokedeath <- mapply(poisson_deviance, m_gof_nm1$n_Strokedeath, m_gof_nm1$n_Strokedeath_fit)

  deviance_all <- m_gof_nm1$deviance_n_death + m_gof_nm1$deviance_n_cvdeath + m_gof_nm1$deviance_n_noncvdeath +
    m_gof_nm1$deviance_n_mi + m_gof_nm1$deviance_n_stroke + m_gof_nm1$deviance_n_tia + m_gof_nm1$deviance_n_hf + m_gof_nm1$deviance_n_strokedeath
  return(deviance_all)
}

# Creating vectors for lower bound and upper bound for parameters
lowb_nm <- c(lowbound_nm["p_HTN_HF"], lowbound_nm["p_MI_HF"], lowbound_nm["p_HF_HF"], lowbound_nm["p_Stroke_HF"], lowbound_nm["p_TIA_HF"], lowbound_nm["p_Othercvd_HF"], lowbound_nm["p_Chroniccvd_HF"], hr_stroke = l_hrparams$hr_stroke * m_lhs_param_values1["hr_multiplier"] * lowbound_nm["hr_stroke_multiplier"])

upperb_nm <- c(upperbound_nm["p_HTN_HF"], upperbound_nm["p_MI_HF"], upperbound_nm["p_HF_HF"], upperbound_nm["p_Stroke_HF"], upperbound_nm["p_TIA_HF"], upperbound_nm["p_Othercvd_HF"], upperbound_nm["p_Chroniccvd_HF"], hr_stroke = l_hrparams$hr_stroke * m_lhs_param_values1["hr_multiplier"] * upperbound_nm["hr_stroke_multiplier"])


numclusters_nm <- makeCluster(10)

# Moving the data frames and objects to all clusters
clusterExport(numclusters_nm, c(
  "v_params_init2_nm", "microsim_nm", "prob_func_nm", "lowb_nm", "upperb_nm",
  "m_lhs_param_values4", "m_gof_nm", "crn", "v.age.monthsfrom60", "m.prob_nm",
  "a_transition_matrix_nm", "background_mortality", "a_state_membership_nm", "n_cycles", "n.i", "cond1.index",
  "cond1"
))

fit_nm_parallel <- parLapply(numclusters_nm, 1:10, function(i) {
  dfoptim::nmkb(par = v_params_init2_nm[i, ], fn = microsim_nm, lower = lowb_nm, upper = upperb_nm, control = list(trace = TRUE))
})

stopCluster(numclusters_nm)

# Saving the output
save(fit_nm_parallel, file = here("Outputs", paste0("fit_nm_parallel_", Sys.Date(), ".RData")))

## Loading the R object ##
##load(here("Outputs", "fit_nm_parallel_2025-05-28.RData"))
## ##

# Identifying the best-fitting values from the Nelder-Mead algorithm
index_fit_nm_parallel <- which.min(sapply(fit_nm_parallel, function(x) x$value))

fit_nm_parallel1 <- fit_nm_parallel[[index_fit_nm_parallel]]

# Extracting the parameter values
fit_nm_best <- fit_nm_parallel1$par

# Combining the two vectors to include all probabilities
v_allparams <- c(as.vector(m_lhs_param_values4), fit_nm_best)

# Adding names
names(v_allparams) <- c(
  "back_multiplier", "p_HTN_MI", "p_HTN_Stroke", "p_HTN_TIA", "p_HTN_Othercvd", "p_MI_MI", "p_MI_Stroke", "p_MI_TIA", "p_MI_Othercvd", "p_HF_MI", "p_HF_Stroke",
  "p_HF_TIA", "p_HF_Othercvd", "p_Stroke_MI", "p_Stroke_Stroke", "p_Stroke_TIA", "p_Stroke_Othercvd", "p_TIA_MI", "p_TIA_Stroke", "p_TIA_TIA",
  "p_TIA_Othercvd", "p_Othercvd_MI", "p_Othercvd_Stroke", "p_Othercvd_TIA", "p_Othercvd_Othercvd", "p_Chroniccvd_MI", "p_Chroniccvd_Stroke", "p_Chroniccvd_TIA",
  "p_HTN_Renal", "p_Renal_Advrenal", "hr_mi", "hr_hf", "hr_tia", "hr_othercvd", "hr_chroniccvd", "hr_renal", "hr_advrenal",
  "p_HTN_HF", "p_MI_HF", "p_HF_HF", "p_Stroke_HF", "p_TIA_HF", "p_Othercvd_HF", "p_Chroniccvd_HF", "hr_stroke"
)


# Converting it to a matrix
m_allparams_nmparallel <- matrix(v_allparams, nrow = 1)

# Adding column names
colnames(m_allparams_nmparallel) <- c(
  "back_multiplier", "p_HTN_MI", "p_HTN_Stroke", "p_HTN_TIA", "p_HTN_Othercvd", "p_MI_MI", "p_MI_Stroke", "p_MI_TIA", "p_MI_Othercvd", "p_HF_MI", "p_HF_Stroke",
  "p_HF_TIA", "p_HF_Othercvd", "p_Stroke_MI", "p_Stroke_Stroke", "p_Stroke_TIA", "p_Stroke_Othercvd", "p_TIA_MI", "p_TIA_Stroke", "p_TIA_TIA",
  "p_TIA_Othercvd", "p_Othercvd_MI", "p_Othercvd_Stroke", "p_Othercvd_TIA", "p_Othercvd_Othercvd", "p_Chroniccvd_MI", "p_Chroniccvd_Stroke", "p_Chroniccvd_TIA",
  "p_HTN_Renal", "p_Renal_Advrenal", "hr_mi", "hr_hf", "hr_tia", "hr_othercvd", "hr_chroniccvd", "hr_renal", "hr_advrenal",
  "p_HTN_HF", "p_MI_HF", "p_HF_HF", "p_Stroke_HF", "p_TIA_HF", "p_Othercvd_HF", "p_Chroniccvd_HF", "hr_stroke"
)


# Saving this matrix of final input parameters
save(m_allparams_nmparallel, file = here("Outputs", paste0("m_allparams_nmparallel_", Sys.Date(), ".RData")))

## Loading the calibrated input parameters##
## load(here("Outputs", "m_allparams_nmparallel_2025-05-28.RData"))
##

# Rerunning the model using the best-fit input parameters####


# Creating an empty transition matrix for rerun of the model using best fit parameters
a_transition_matrix_bestfit <- array(NA_real_, dim = c(n_states, n_states, n_cycles), dimnames = list(
  from = cond1,
  to = cond1,
  cycle = 1:n_cycles
))

# Creating an empty array for the simulation cohort
a_state_membership_bestfit <- array(NA_real_,
  dim = c(n_cycles, n_states),
  dimnames = list(
    cycle = 1:n_cycles,
    state = cond1
  )
)

# Creating a matrix for states
m.prob_bestfit <- matrix(
  nrow = n.i, ncol = n_cycles,
  dimnames = list(
    paste("ind", 1:n.i, sep = " "),
    paste("cycle", 1:n_cycles, sep = " ")
  )
)

# Initializing the first cycle as HTN state
m.prob_bestfit[, 1] <- "HTN"

# Goodness of fit matrix
m_gof_bestfit <- matrix(nrow = 1, ncol = n_lhs_targets * 2)

colnames(m_gof_bestfit)[9:16] <- paste0(v_lhs_targets, "_fit")

colnames(m_gof_bestfit)[1:8] <- v_lhs_targets

# Values for calibration targets
m_gof_bestfit[, "n_Death"] <- 45000 # 169 cycles
m_gof_bestfit[, "n_Cvdeath"] <- 24300 # 168 cycles
m_gof_bestfit[, "n_Noncvdeath"] <- 20700 # 168 cycles
m_gof_bestfit[, "n_Strokedeath"] <- 5000 # 264 cycles
m_gof_bestfit[, "n_MI"] <- 5054 # 60 cycles
m_gof_bestfit[, "n_HF"] <- 5483 # 60 cycles
m_gof_bestfit[, "n_Stroke"] <- 8200 # 60 cycles
m_gof_bestfit[, "n_TIA"] <- 4229 # 60 cycles

# Function for model and goodness of fit
prob_func_bestfit <- function(params) {
  ## Filling the transition probability values##

  with(as.list(as.data.frame(params)), {
    # From HTN to Death across all cycles
    # HTN to non-cardiovascular death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HTN", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HTN", "Noncvdeath", x] > 1] <- 1
    }

    # From HTN to Other states
    a_transition_matrix_bestfit["HTN", "MI", ] <- p_HTN_MI
    a_transition_matrix_bestfit["HTN", "HF", ] <- p_HTN_HF
    a_transition_matrix_bestfit["HTN", "Stroke", ] <- p_HTN_Stroke
    a_transition_matrix_bestfit["HTN", "TIA", ] <- p_HTN_TIA
    a_transition_matrix_bestfit["HTN", "Renal", ] <- p_HTN_Renal
    a_transition_matrix_bestfit["HTN", "Othercvd", ] <- p_HTN_Othercvd
    a_transition_matrix_bestfit["HTN", "HTN", ] <- 1 - apply(a_transition_matrix_bestfit["HTN", , ], 2, sum, na.rm = TRUE)

    # From MI to all other states
    # MI to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI", "Noncvdeath", x] > 1] <- 1
    }

    # MI to MI2
    a_transition_matrix_bestfit["MI", "MI2", ] <- 1 - apply(a_transition_matrix_bestfit["MI", , ], 2, sum, na.rm = TRUE)

    # From MI2 to Other states
    # MI2 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI2", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI2 to Non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI2", "Noncvdeath", x] > 1] <- 1
    }

    # MI2 to MI3
    a_transition_matrix_bestfit["MI2", "MI3", ] <- 1 - apply(a_transition_matrix_bestfit["MI2", , ], 2, sum, na.rm = TRUE)

    # From MI3 to Other states
    # MI3 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI3", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI3 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI3", "Noncvdeath", x] > 1] <- 1
    }

    # MI3 to MI4
    a_transition_matrix_bestfit["MI3", "MI4", ] <- 1 - apply(a_transition_matrix_bestfit["MI3", , ], 2, sum, na.rm = TRUE)

    # From MI4 to Other states
    # MI4 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI4", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI4 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI4", "Noncvdeath", x] > 1] <- 1
    }

    # MI4 to MI5
    a_transition_matrix_bestfit["MI4", "MI5", ] <- 1 - apply(a_transition_matrix_bestfit["MI4", , ], 2, sum, na.rm = TRUE)

    # From MI5 to Other states
    # MI5 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI5", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI5 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI5", "Noncvdeath", x] > 1] <- 1
    }

    # MI5 to MI6
    a_transition_matrix_bestfit["MI5", "MI6", ] <- 1 - apply(a_transition_matrix_bestfit["MI5", , ], 2, sum, na.rm = TRUE)

    # From MI6 to Other states
    # MI6 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI6", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI6 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI6", "Noncvdeath", x] > 1] <- 1
    }

    # MI6 to MI7
    a_transition_matrix_bestfit["MI6", "MI7", ] <- 1 - apply(a_transition_matrix_bestfit["MI6", , ], 2, sum, na.rm = TRUE)

    # From MI7 to Other states
    # MI7 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI7", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI7 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI7", "Noncvdeath", x] > 1] <- 1
    }

    # MI7 to MI8
    a_transition_matrix_bestfit["MI7", "MI8", ] <- 1 - apply(a_transition_matrix_bestfit["MI7", , ], 2, sum, na.rm = TRUE)

    # From MI8 to Other states
    # MI8 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI8", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI8 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI8", "Noncvdeath", x] > 1] <- 1
    }

    # MI8 to MI9
    a_transition_matrix_bestfit["MI8", "MI9", ] <- 1 - apply(a_transition_matrix_bestfit["MI8", , ], 2, sum, na.rm = TRUE)

    # From MI9 to Other states
    # MI9 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI9", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI9 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI9", "Noncvdeath", x] > 1] <- 1
    }

    # MI9 to MI10
    a_transition_matrix_bestfit["MI9", "MI10", ] <- 1 - apply(a_transition_matrix_bestfit["MI9", , ], 2, sum, na.rm = TRUE)

    # From MI10 to Other states
    # MI10 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI10", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI10 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI10", "Noncvdeath", x] > 1] <- 1
    }

    # MI10 to MI11
    a_transition_matrix_bestfit["MI10", "MI11", ] <- 1 - apply(a_transition_matrix_bestfit["MI10", , ], 2, sum, na.rm = TRUE)

    # From MI11 to Other states
    # MI11 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI11", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # MI11 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI11", "Noncvdeath", x] > 1] <- 1
    }

    # MI11 to MI12
    a_transition_matrix_bestfit["MI11", "MI12", ] <- 1 - apply(a_transition_matrix_bestfit["MI11", , ], 2, sum, na.rm = TRUE)

    # From MI12 to Other states
    # MI12 to MI death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI12", "MIdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_mi))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI12", "MIdeath", x] > 1] <- 1
    }

    # MI12 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["MI12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["MI12", "Noncvdeath", x] > 1] <- 1
    }

    a_transition_matrix_bestfit["MI12", "MI", ] <- p_MI_MI
    a_transition_matrix_bestfit["MI12", "HF", ] <- p_MI_HF
    a_transition_matrix_bestfit["MI12", "Stroke", ] <- p_MI_Stroke
    a_transition_matrix_bestfit["MI12", "TIA", ] <- p_MI_TIA
    a_transition_matrix_bestfit["MI12", "Othercvd", ] <- p_MI_Othercvd
    a_transition_matrix_bestfit["MI12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_bestfit["MI12", , ], 2, sum, na.rm = TRUE)

    # From HF to Other states
    # HF to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF", "HFdeath", x] > 1] <- 1
    }

    # HF to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF", "Noncvdeath", x] > 1] <- 1
    }

    # HF to HF2
    a_transition_matrix_bestfit["HF", "HF2", ] <- 1 - apply(a_transition_matrix_bestfit["HF", , ], 2, sum, na.rm = TRUE)

    # From HF2 to Other states
    # HF2 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF2", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF2", "HFdeath", x] > 1] <- 1
    }

    # HF2 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF2", "Noncvdeath", x] > 1] <- 1
    }

    # HF2 to HF3
    a_transition_matrix_bestfit["HF2", "HF3", ] <- 1 - apply(a_transition_matrix_bestfit["HF2", , ], 2, sum, na.rm = TRUE)

    # From HF3 to Other states
    # HF3 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF3", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF3", "HFdeath", x] > 1] <- 1
    }

    # HF3 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF3", "Noncvdeath", x] > 1] <- 1
    }

    # HF3 to HF4
    a_transition_matrix_bestfit["HF3", "HF4", ] <- 1 - apply(a_transition_matrix_bestfit["HF3", , ], 2, sum, na.rm = TRUE)

    # From HF4 to Other states
    # HF4 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF4", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF4", "HFdeath", x] > 1] <- 1
    }

    # HF4 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF4", "Noncvdeath", x] > 1] <- 1
    }

    # HF4 to HF5
    a_transition_matrix_bestfit["HF4", "HF5", ] <- 1 - apply(a_transition_matrix_bestfit["HF4", , ], 2, sum, na.rm = TRUE)

    # From HF5 to Other states
    # HF5 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF5", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF5", "HFdeath", x] > 1] <- 1
    }

    # HF5 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF5", "Noncvdeath", x] > 1] <- 1
    }

    # HF5 to HF6
    a_transition_matrix_bestfit["HF5", "HF6", ] <- 1 - apply(a_transition_matrix_bestfit["HF5", , ], 2, sum, na.rm = TRUE)

    # From HF6 to Other states
    # HF6 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF6", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF6", "HFdeath", x] > 1] <- 1
    }

    # HF6 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF6", "Noncvdeath", x] > 1] <- 1
    }

    # HF6 to HF7
    a_transition_matrix_bestfit["HF6", "HF7", ] <- 1 - apply(a_transition_matrix_bestfit["HF6", , ], 2, sum, na.rm = TRUE)

    # From HF7 to Other states
    # HF7 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF7", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF7", "HFdeath", x] > 1] <- 1
    }

    # HF7 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF7", "Noncvdeath", x] > 1] <- 1
    }

    # HF7 to HF8
    a_transition_matrix_bestfit["HF7", "HF8", ] <- 1 - apply(a_transition_matrix_bestfit["HF7", , ], 2, sum, na.rm = TRUE)

    # From HF8 to Other states
    # HF8 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF8", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF8", "HFdeath", x] > 1] <- 1
    }

    # HF8 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF8", "Noncvdeath", x] > 1] <- 1
    }

    # HF8 to HF9
    a_transition_matrix_bestfit["HF8", "HF9", ] <- 1 - apply(a_transition_matrix_bestfit["HF8", , ], 2, sum, na.rm = TRUE)

    # From HF9 to Other states
    # HF9 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF9", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF9", "HFdeath", x] > 1] <- 1
    }

    # HF9 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF9", "Noncvdeath", x] > 1] <- 1
    }

    # HF9 to HF10
    a_transition_matrix_bestfit["HF9", "HF10", ] <- 1 - apply(a_transition_matrix_bestfit["HF9", , ], 2, sum, na.rm = TRUE)

    # From HF10 to Other states
    # HF10 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF10", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF10", "HFdeath", x] > 1] <- 1
    }

    # HF10 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF10", "Noncvdeath", x] > 1] <- 1
    }

    # HF10 to HF11
    a_transition_matrix_bestfit["HF10", "HF11", ] <- 1 - apply(a_transition_matrix_bestfit["HF10", , ], 2, sum, na.rm = TRUE)

    # From HF11 to Other states
    # HF11 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF11", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF11", "HFdeath", x] > 1] <- 1
    }

    # HF11 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF11", "Noncvdeath", x] > 1] <- 1
    }

    # HF11 to HF12
    a_transition_matrix_bestfit["HF11", "HF12", ] <- 1 - apply(a_transition_matrix_bestfit["HF11", , ], 2, sum, na.rm = TRUE)

    # From HF12 to Other states
    # HF12 to HF death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF12", "HFdeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_hf))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF12", "HFdeath", x] > 1] <- 1
    }

    # HF12 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["HF12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["HF12", "Noncvdeath", x] > 1] <- 1
    }

    a_transition_matrix_bestfit["HF12", "MI", ] <- p_HF_MI
    a_transition_matrix_bestfit["HF12", "HF", ] <- p_HF_HF
    a_transition_matrix_bestfit["HF12", "Stroke", ] <- p_HF_Stroke
    a_transition_matrix_bestfit["HF12", "TIA", ] <- p_HF_TIA
    a_transition_matrix_bestfit["HF12", "Othercvd", ] <- p_HF_Othercvd
    a_transition_matrix_bestfit["HF12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_bestfit["HF12", , ], 2, sum, na.rm = TRUE)

    # From Stroke to Other states
    # Stroke to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke to Stroke2
    a_transition_matrix_bestfit["Stroke", "Stroke2", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke", , ], 2, sum, na.rm = TRUE)

    # From Stroke2 to Other states
    # Stroke2 to Non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke2", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke2 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke2", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke2 to stroke3
    a_transition_matrix_bestfit["Stroke2", "Stroke3", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke2", , ], 2, sum, na.rm = TRUE)

    # From Stroke3 to Other states
    # Stroke3 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke3", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke3 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke3", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke3 to stroke4
    a_transition_matrix_bestfit["Stroke3", "Stroke4", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke3", , ], 2, sum, na.rm = TRUE)

    # From Stroke4 to Other states
    # Stroke4 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke4", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke4 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke4", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke4 to Stroke5
    a_transition_matrix_bestfit["Stroke4", "Stroke5", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke4", , ], 2, sum, na.rm = TRUE)

    # From Stroke5 to Other states
    # Stroke5 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke5", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke5 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke5", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke5 to Stroke6
    a_transition_matrix_bestfit["Stroke5", "Stroke6", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke5", , ], 2, sum, na.rm = TRUE)

    # From Stroke6 to Other states
    # Stroke6 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke6", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke6 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke6", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke6 to Stroke7
    a_transition_matrix_bestfit["Stroke6", "Stroke7", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke6", , ], 2, sum, na.rm = TRUE)

    # From Stroke7 to Other states

    # Stroke7 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke7", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke7 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke7", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke7 to Stroke8
    a_transition_matrix_bestfit["Stroke7", "Stroke8", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke7", , ], 2, sum, na.rm = TRUE)

    # From Stroke8 to Other states
    # Stroke8 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke8", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke8 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke8", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke8 to Stroke9
    a_transition_matrix_bestfit["Stroke8", "Stroke9", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke8", , ], 2, sum, na.rm = TRUE)

    # From Stroke9 to Other states
    # Stroke9 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke9", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke9 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke9", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke9 to Stroke10
    a_transition_matrix_bestfit["Stroke9", "Stroke10", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke9", , ], 2, sum, na.rm = TRUE)

    # From Stroke10 to Other states

    # Stroke10 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke10", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke10 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke10", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke10 to Stroke11
    a_transition_matrix_bestfit["Stroke10", "Stroke11", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke10", , ], 2, sum, na.rm = TRUE)

    # From Stroke11 to Other states
    # Stroke11 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke11", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke11 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke11", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # Stroke11 to Stroke12
    a_transition_matrix_bestfit["Stroke11", "Stroke12", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke11", , ], 2, sum, na.rm = TRUE)

    # From Stroke12 to Other states
    # Stroke12 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Stroke12", "Noncvdeath", x] > 1] <- 1
    }

    # Stroke12 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Stroke12", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_stroke))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    a_transition_matrix_bestfit["Stroke12", "MI", ] <- p_Stroke_MI
    a_transition_matrix_bestfit["Stroke12", "HF", ] <- p_Stroke_HF
    a_transition_matrix_bestfit["Stroke12", "Stroke", ] <- p_Stroke_Stroke
    a_transition_matrix_bestfit["Stroke12", "TIA", ] <- p_Stroke_TIA
    a_transition_matrix_bestfit["Stroke12", "Othercvd", ] <- p_Stroke_Othercvd
    a_transition_matrix_bestfit["Stroke12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_bestfit["Stroke12", , ], 2, sum, na.rm = TRUE)

    # From TIA to Other states
    # TIA to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA", "Noncvdeath", x] > 1] <- 1
    }

    # TIA to TIA2
    a_transition_matrix_bestfit["TIA", "TIA2", ] <- 1 - apply(a_transition_matrix_bestfit["TIA", , ], 2, sum, na.rm = TRUE)

    # From TIA2 to Other states
    # TIA2 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA2", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA2 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA2", "Noncvdeath", x] > 1] <- 1
    }

    # TIA2 to TIA3
    a_transition_matrix_bestfit["TIA2", "TIA3", ] <- 1 - apply(a_transition_matrix_bestfit["TIA2", , ], 2, sum, na.rm = TRUE)

    # From TIA3 to Other states
    # TIA3 to Stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA3", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA3 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA3", "Noncvdeath", x] > 1] <- 1
    }

    # TIA3 to TIa4
    a_transition_matrix_bestfit["TIA3", "TIA4", ] <- 1 - apply(a_transition_matrix_bestfit["TIA3", , ], 2, sum, na.rm = TRUE)

    # From TIA4 to Other states
    # TIA4 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA4", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA4 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA4", "Noncvdeath", x] > 1] <- 1
    }

    # TIA4 to TIA5
    a_transition_matrix_bestfit["TIA4", "TIA5", ] <- 1 - apply(a_transition_matrix_bestfit["TIA4", , ], 2, sum, na.rm = TRUE)

    # From TIA5 to Other states
    # TIA5 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA5", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA5 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA5", "Noncvdeath", x] > 1] <- 1
    }

    # TIA5 to TIA6
    a_transition_matrix_bestfit["TIA5", "TIA6", ] <- 1 - apply(a_transition_matrix_bestfit["TIA5", , ], 2, sum, na.rm = TRUE)

    # From TIA6 to Other states
    # TIA6 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA6", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA6 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA6", "Noncvdeath", x] > 1] <- 1
    }

    # TIA6 to TIA7
    a_transition_matrix_bestfit["TIA6", "TIA7", ] <- 1 - apply(a_transition_matrix_bestfit["TIA6", , ], 2, sum, na.rm = TRUE)

    # From TIA7 to Other states
    # TIA7 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA7", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA7  to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA7", "Noncvdeath", x] > 1] <- 1
    }

    # TIA7 to TIA8
    a_transition_matrix_bestfit["TIA7", "TIA8", ] <- 1 - apply(a_transition_matrix_bestfit["TIA7", , ], 2, sum, na.rm = TRUE)

    # From TIA8 to Other states
    # TIA8 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA8", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA8 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA8", "Noncvdeath", x] > 1] <- 1
    }

    # TIA8 to TIA9
    a_transition_matrix_bestfit["TIA8", "TIA9", ] <- 1 - apply(a_transition_matrix_bestfit["TIA8", , ], 2, sum, na.rm = TRUE)

    # From TIA9 to Other states
    # TIA9 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA9", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA9 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA9", "Noncvdeath", x] > 1] <- 1
    }

    # TIA9 to TIA10
    a_transition_matrix_bestfit["TIA9", "TIA10", ] <- 1 - apply(a_transition_matrix_bestfit["TIA9", , ], 2, sum, na.rm = TRUE)

    # From TIA10 to Other states
    # TIA10 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA10", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA10 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA10", "Noncvdeath", x] > 1] <- 1
    }

    # TIA10 to TIA11
    a_transition_matrix_bestfit["TIA10", "TIA11", ] <- 1 - apply(a_transition_matrix_bestfit["TIA10", , ], 2, sum, na.rm = TRUE)

    # From TIA11 to Other states
    # TIA11 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA11", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA11 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA11", "Noncvdeath", x] > 1] <- 1
    }

    # TIA11 to TIA12
    a_transition_matrix_bestfit["TIA11", "TIA12", ] <- 1 - apply(a_transition_matrix_bestfit["TIA11", , ], 2, sum, na.rm = TRUE)

    # From TIA12 to Other states
    # TIA12 to stroke death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA12", "Strokedeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_tia))) - (1 - exp(-(background_mortality$mnth_rate[x])))
    }

    # TIA12 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["TIA12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA12", "Noncvdeath", x] > 1] <- 1
    }

    a_transition_matrix_bestfit["TIA12", "MI", ] <- p_TIA_MI
    a_transition_matrix_bestfit["TIA12", "HF", ] <- p_TIA_HF
    a_transition_matrix_bestfit["TIA12", "Stroke", ] <- p_TIA_Stroke
    a_transition_matrix_bestfit["TIA12", "TIA", ] <- p_TIA_TIA
    a_transition_matrix_bestfit["TIA12", "Othercvd", ] <- p_TIA_Othercvd
    a_transition_matrix_bestfit["TIA12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_bestfit["TIA12", , ], 2, sum, na.rm = TRUE)

    # From Othercvd to other states
    # Othercvd to non-cv death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd", "Othercvddeath", x] > 1] <- 1
    }
    # Othercvd to Othercvd2
    a_transition_matrix_bestfit["Othercvd", "Othercvd2", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd", , ], 2, sum, na.rm = TRUE)

    # From Othercvd2 to Other states
    # Othercvd2 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd2", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd2", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd2 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd2", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd2", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd2 to Othercvd3
    a_transition_matrix_bestfit["Othercvd2", "Othercvd3", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd2", , ], 2, sum, na.rm = TRUE)

    # From Othercvd3 to Other states
    # Othercvd3 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd3", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd3", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd3 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd3", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd3", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd3 to Othercvd4
    a_transition_matrix_bestfit["Othercvd3", "Othercvd4", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd3", , ], 2, sum, na.rm = TRUE)

    # From Othercvd4 to Other states
    # Othercvd4 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd4", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd4", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd4 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd4", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd4", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd4 to Othercvd5
    a_transition_matrix_bestfit["Othercvd4", "Othercvd5", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd4", , ], 2, sum, na.rm = TRUE)

    # From Othercvd5 to Other states
    # Othercvd5 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd5", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd5", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd5 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd5", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd5", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd5 to Othercvd6
    a_transition_matrix_bestfit["Othercvd5", "Othercvd6", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd5", , ], 2, sum, na.rm = TRUE)

    # From Othercvd6 to Other states
    # Othercvd6 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd6", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd6", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd6 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd6", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd6", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd6 to Othercvd7
    a_transition_matrix_bestfit["Othercvd6", "Othercvd7", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd6", , ], 2, sum, na.rm = TRUE)

    # From Othercvd7 to Other states
    # Othercvd7 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd7", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd7", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd7 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd7", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd7", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd7 to Othercvd8
    a_transition_matrix_bestfit["Othercvd7", "Othercvd8", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd7", , ], 2, sum, na.rm = TRUE)

    # From Othercvd8 to Other states
    # Othercvd8 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd8", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd8", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd8 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd8", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd8", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd8 to Othercvd9
    a_transition_matrix_bestfit["Othercvd8", "Othercvd9", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd8", , ], 2, sum, na.rm = TRUE)

    # From Othercvd9 to Other states
    # Othercvd9 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd9", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd9", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd9 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd9", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd9", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd9 to Othercvd10
    a_transition_matrix_bestfit["Othercvd9", "Othercvd10", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd9", , ], 2, sum, na.rm = TRUE)

    # From Othercvd10 to Other states
    # Othercvd10 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd10", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["TIA10", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd10 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd10", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd10", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd10 to Othercvd11
    a_transition_matrix_bestfit["Othercvd10", "Othercvd11", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd10", , ], 2, sum, na.rm = TRUE)

    # From Othercvd11 to Other states
    # Othercvd11 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd11", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd11", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd11 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd11", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd11", "Noncvdeath", x] > 1] <- 1
    }

    # Othercvd11 to Othercvd12
    a_transition_matrix_bestfit["Othercvd11", "Othercvd12", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd11", , ], 2, sum, na.rm = TRUE)

    # From Othercvd12 to Other states
    # Othercvd12 to Othercvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd12", "Othercvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_othercvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd12", "Othercvddeath", x] > 1] <- 1
    }

    # Othercvd12 to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Othercvd12", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Othercvd12", "Noncvdeath", x] > 1] <- 1
    }

    a_transition_matrix_bestfit["Othercvd12", "MI", ] <- p_Othercvd_MI
    a_transition_matrix_bestfit["Othercvd12", "HF", ] <- p_Othercvd_HF
    a_transition_matrix_bestfit["Othercvd12", "Stroke", ] <- p_Othercvd_Stroke
    a_transition_matrix_bestfit["Othercvd12", "TIA", ] <- p_Othercvd_TIA
    a_transition_matrix_bestfit["Othercvd12", "Othercvd", ] <- p_Othercvd_Othercvd
    a_transition_matrix_bestfit["Othercvd12", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_bestfit["Othercvd12", , ], 2, sum, na.rm = TRUE)


    # From Renal disease to Other states
    # Renal to non-CV death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Renal", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier * hr_renal
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Renal", "Noncvdeath", x] > 1] <- 1
    }

    a_transition_matrix_bestfit["Renal", "Advrenal", ] <- p_Renal_Advrenal
    a_transition_matrix_bestfit["Renal", "Renal", ] <- 1 - apply(a_transition_matrix_bestfit["Renal", , ], 2, sum, na.rm = TRUE)

    # From Chroniccvd to Other states
    # Chronic cvd to non-cv death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Chroniccvd", "Noncvdeath", x] <- background_mortality$mnth_prob[x] * back_multiplier
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Chroniccvd", "Noncvdeath", x] > 1] <- 1
    }

    # Chronic cvd to Chroniccvd death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Chroniccvd", "Chroniccvddeath", x] <- (1 - exp(-(background_mortality$mnth_rate[x] * hr_chroniccvd))) - (1 - exp(-(background_mortality$mnth_rate[x])))
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Chroniccvd", "Chroniccvddeath", x] > 1] <- 1
    }

    a_transition_matrix_bestfit["Chroniccvd", "MI", ] <- p_Chroniccvd_MI
    a_transition_matrix_bestfit["Chroniccvd", "HF", ] <- p_Chroniccvd_HF
    a_transition_matrix_bestfit["Chroniccvd", "Stroke", ] <- p_Chroniccvd_Stroke
    a_transition_matrix_bestfit["Chroniccvd", "TIA", ] <- p_Chroniccvd_TIA
    a_transition_matrix_bestfit["Chroniccvd", "Chroniccvd", ] <- 1 - apply(a_transition_matrix_bestfit["Chroniccvd", , ], 2, sum, na.rm = TRUE)

    # From Advanced renal to Other states

    # Advanced renal to non-cv death
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit["Advrenal", "Noncvdeath", x] <- (background_mortality$mnth_prob[x] * back_multiplier * hr_advrenal)
      a_transition_matrix_bestfit[a_transition_matrix_bestfit["Advrenal", "Noncvdeath", x] > 1] <- 1
    }
    a_transition_matrix_bestfit["Advrenal", "Advrenal", ] <- 1 - apply(a_transition_matrix_bestfit["Advrenal", , ], 2, sum, na.rm = TRUE)

    # From MI death to Death
    a_transition_matrix_bestfit["MIdeath", "Death", ] <- 1

    # From HF death to Death
    a_transition_matrix_bestfit["HFdeath", "Death", ] <- 1

    # From stroke death to Death
    a_transition_matrix_bestfit["Strokedeath", "Death", ] <- 1

    # From Othercvd death to Death
    a_transition_matrix_bestfit["Othercvddeath", "Death", ] <- 1

    # From non-CV death to Death
    a_transition_matrix_bestfit["Noncvdeath", "Death", ] <- 1

    # From chroniccvd death to Death
    a_transition_matrix_bestfit["Chroniccvddeath", "Death", ] <- 1

    # From Death to Other states
    a_transition_matrix_bestfit["Death", "Death", ] <- 1

    # Replace NAs with 0 in the transition_matrix
    a_transition_matrix_bestfit[is.na(a_transition_matrix_bestfit)] <- 0

    # Standardizing the probabilities so that they sum up to 1
    for (x in 1:n_cycles) {
      a_transition_matrix_bestfit[, , x] <- a_transition_matrix_bestfit[, , x] / (apply(a_transition_matrix_bestfit[, , x], MARGIN = 1, FUN = sum, na.rm = TRUE))
    }
    return(a_transition_matrix_bestfit)
  })
}

# Creating the function for microsimulation model with common random numbers####
microsim_bestfit <- function(params) {
  # Calling the function for transition probability matrix
  m_transitionprobs <- prob_func_bestfit(params)

  for (n in 1:n.i) {
    set.seed(crn[n])
    v.age.monthsstart <- v.age.monthsfrom60[n]
    a_state_membership_bestfit[v.age.monthsstart, ] <- c(
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
    for (i in (v.age.monthsstart + 1):n_cycles) {
      a_state_membership_bestfit[i, ] <- a_state_membership_bestfit[i - 1, ] %*% m_transitionprobs[, , i - 1]
      index.state <- sample(cond1.index, prob = a_state_membership_bestfit[i, ], size = 1)
      m.prob_bestfit[n, i - v.age.monthsstart + 1] <- cond1[index.state]
      a_state_membership_bestfit[i, ] <- c(
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      )
      a_state_membership_bestfit[i, index.state] <- 1
    }
  }

  # Replacing all NAs with Death
  m.prob1_bestfit <- m.prob_bestfit
  m.prob1_bestfit[is.na(m.prob1_bestfit)] <- "Death"

  # Converting to dataframe
  m.prob2_bestfit <- as.data.frame(m.prob1_bestfit)

  # Calculating how many patients died in cycles 1 - 169
  m_gof_bestfit[, "n_Death_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:169)] == "Death", 1, any), ])

  # Calculating how many patients had CV death in cycles 1 - 168
  m_gof_bestfit[, "n_Cvdeath_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:168)] == "MIdeath", 1, any), ]) +
    nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:168)] == "HFdeath", 1, any), ]) +
    nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:168)] == "Strokedeath", 1, any), ]) +
    nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:168)] == "Othercvddeath", 1, any), ]) +
    nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:168)] == "Chroniccvddeath", 1, any), ])

  # Calculating how many patients had non-CV death in cycles 1 - 168
  m_gof_bestfit[, "n_Noncvdeath_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:168)] == "Noncvdeath", 1, any), ])

  # Calculating how many patients had stroke death in cycles 1 - 264
  m_gof_bestfit[, "n_Strokedeath_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:264)] == "Strokedeath", 1, any), ])

  # Calculating how many patients had stroke in cycles 1- 60
  m_gof_bestfit[, "n_Stroke_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:60)] == "Stroke", 1, any), ])

  # Calculating how many patients had TIA in cycles 1 - 60
  m_gof_bestfit[, "n_TIA_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:60)] == "TIA", 1, any), ])

  # Calculating how many patients had MI in cycles 1 - 60
  m_gof_bestfit[, "n_MI_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:60)] == "MI", 1, any), ])

  # Calculating how many patients had HF in cycles 1 - 60
  m_gof_bestfit[, "n_HF_fit"] <- nrow(m.prob2_bestfit[apply(m.prob2_bestfit[c(1:60)] == "HF", 1, any), ])

  return(list(m.prob2_bestfit, m_gof_bestfit))
}

# Output the function
bestfit_results_nm <- microsim_bestfit(m_allparams_nmparallel)

# Saving the R object
save(bestfit_results_nm, file = here("Outputs", paste0("bestfit_results_nm_", Sys.Date(), ".RData")))

## Loading the R object##
## load(here("Outputs", "bestfit_results_nm_2025-06-01.RData"))
##

# Extracting the dataframe with number of patients with different outcomes used in calibration
m_gof_bestfit <- bestfit_results_nm[[2]]

# Extracting the transition matrix across cycles
matrix_allcycles <- (bestfit_results_nm[[1]])

# Creating a data frame for cumulative incidence rates by year for 10 years for MI, CHF, Stroke and TIA
df_for_plot <- data.frame(
  MI = c(0, nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "MI", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "MI", 1, any), ])),
  CHF = c(0, nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "HF", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "HF", 1, any), ])),
  Stroke = c(0, nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "Stroke", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "Stroke", 1, any), ])),
  TIA = c(0, nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "TIA", 1, any), ]), nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "TIA", 1, any), ])),
  time = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
)

# Transforming data from wide to long
df_for_plot1 <- df_for_plot %>% pivot_longer(cols = c("MI", "CHF", "Stroke", "TIA"), names_to = "Type", values_to = "Inc.rates")

# Plotting the incidence rates
inc_rates_plot <- ggplot(data = df_for_plot1, aes(x = time, y = Inc.rates, color = Type)) +
  geom_line() +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  scale_color_viridis(discrete = TRUE) +
  theme_classic() +
  xlab("Time in years") +
  ylab("Cumulative Incidence per 100,000")

# Saving the plot (as .svg)
ggsave(
  filename = here("Outputs", paste0("inc_rates_", Sys.Date(), ".svg")),
  plot = inc_rates_plot, width = 10, height = 10
)

# Saving the plot (as .png)
ggsave(
  filename = here("Outputs", paste0("inc_rates_", Sys.Date(), ".png")),
  plot = inc_rates_plot, width = 8, height = 8
)

# Looking at observed (digitized) rates and simulated rates of all-cause mortality and
# cardiovascular mortality
# Importing the digitized values
df_mortality_rates <- readxl::read_excel(here("Data", "Extracted points_KM curves.xlsx"))

df_mortality_rates$Time <- round(df_mortality_rates$Time, digits = 0)

# Calculating CV survival probability from the simulations
v_Cvdeaths_year <- c(nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:12)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:24)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:36)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:48)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:60)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:72)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:84)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:96)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:108)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:120)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:132)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:132)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:132)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:132)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:132)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:144)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:144)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:144)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:144)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:144)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:156)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:156)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:156)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:156)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:156)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:168)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:168)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:168)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:168)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:168)] == "Chroniccvddeath", 1, any), ]), 
                     nrow(matrix_allcycles[apply(matrix_allcycles[c(1:180)] == "MIdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:180)] == "HFdeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:180)] == "Strokedeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:180)] == "Othercvddeath", 1, any), ]) + nrow(matrix_allcycles[apply(matrix_allcycles[c(1:180)] == "Chroniccvddeath", 1, any), ]))

# Converting to survival probability
v_Cvdeaths_year1 <- (1 - (v_Cvdeaths_year / 100000))

# Calculating overall survival probability from the simulations
v_alldeaths_year <- c(nrow(matrix_allcycles[apply(matrix_allcycles[c(1:13)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:25)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:37)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:49)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:61)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:73)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:85)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:97)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:109)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:121)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:133)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:145)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:157)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:169)] == "Death", 1, any), ]), 
                      nrow(matrix_allcycles[apply(matrix_allcycles[c(1:181)] == "Death", 1, any), ]))

# Converting to survival probability
v_alldeaths_year1 <- (1 - (v_alldeaths_year / 100000))

# Creating a dataframe from simulated all-cause mortality
df_allcausemortality_rates_s <- data.frame(
  Time = c(365, 365 * 2, 365 * 3, 365 * 4, 365 * 5, 365 * 6, 365 * 7, 365 * 8, 365 * 9, 365 * 10, 365 * 11, 365 * 12, 365 * 13, 365 * 14, 365 * 15),
  Survival = v_alldeaths_year1,
  Type = rep("Overall survival", 15),
  Source = rep("Simulated", 15)
)

# Creating a dataframe from simulated cardiovascular mortality
df_cvmortality_rates_s <- data.frame(
  Time = c(365, 365 * 2, 365 * 3, 365 * 4, 365 * 5, 365 * 6, 365 * 7, 365 * 8, 365 * 9, 365 * 10, 365 * 11, 365 * 12, 365 * 13, 365 * 14, 365 * 15),
  Survival = v_Cvdeaths_year1,
  Type = rep("Cardiovascular survival", 15),
  Source = rep("Simulated", 15)
)

# Combining observed and simulated dataframes
df_mortality_combined <- rbind(df_mortality_rates, df_allcausemortality_rates_s, df_cvmortality_rates_s)

# Plot
mortality_plots <- ggplot(data = df_mortality_combined, aes(x = Time, y = Survival, color = Source)) +
  geom_line() +
  facet_wrap(~Type) +
  scale_color_viridis(discrete = TRUE) +
  theme_classic() +
  xlab("Time in days") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

# Saving the plot as png file
ggsave(
  filename = here("Outputs", paste0("mortality_plots_", Sys.Date(), ".png")),
  plot = mortality_plots, width = 8, height = 8
)

# Saving the plot as svg file
ggsave(
  filename = here("Outputs", paste0("mortality_plots_", Sys.Date(), ".svg")),
  plot = mortality_plots, width = 10, height = 10
)


# External validation of the model####

# Looking at recurrent rates of stroke and MI
matrix_allcycles1 <- matrix_allcycles

# Creating a function for minimum cycle number for MI and death states
firstoccur_func <- function(row, type) {
  min(which(row == type))
}

# Applying the function
matrix_allcycles1$mi.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "MI")
matrix_allcycles1$mi12.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "MI12")
matrix_allcycles1$stroke.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Stroke")
matrix_allcycles1$stroke2.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Stroke2")
matrix_allcycles1$stroke12.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Stroke12")
matrix_allcycles1$hf.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "HF")
matrix_allcycles1$hf12.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "HF12")
matrix_allcycles1$mideath.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "MIdeath")
matrix_allcycles1$strokedeath.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Strokedeath")
matrix_allcycles1$hfdeath.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "HFdeath")
matrix_allcycles1$othercvddeath.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Othercvddeath")
matrix_allcycles1$chroniccvddeath.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Chroniccvddeath")
matrix_allcycles1$noncvdeath.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Noncvdeath")
matrix_allcycles1$death.index <- apply(matrix_allcycles1, 1, firstoccur_func, type = "Death")

# Creating a function for 2nd minimum cycle number for MI and death states
secondoccur_func <- function(row, type) {
  positions <- which(row == type)
  return(positions[2])
}

# Applying the function for MI
matrix_allcycles1$mi.2ndindex <- apply(matrix_allcycles1, 1, secondoccur_func, type = "MI")
matrix_allcycles1$mi12.2ndindex <- apply(matrix_allcycles1, 1, secondoccur_func, type = "MI12")

# Subsetting those individuals that had MI event
matrix_allcycles1_mi <- matrix_allcycles1[!is.infinite(matrix_allcycles1$mi.index), ]

# Creating a minimum death column
matrix_allcycles1_mi$mindeath.index <- pmin(matrix_allcycles1_mi$mideath.index, matrix_allcycles1_mi$strokedeath.index, matrix_allcycles1_mi$hfdeath.index,
                                            matrix_allcycles1_mi$othercvddeath.index, matrix_allcycles1_mi$chroniccvddeath.index, matrix_allcycles1_mi$noncvdeath.index,
                                            matrix_allcycles1_mi$death.index)

# Looking at time to death
# 1 year (12 months)
matrix_allcycles1_mi$death.dummy.12mnth <- ifelse(matrix_allcycles1_mi$mindeath.index <= (matrix_allcycles1_mi$mi.index + 11), 1, 0)
matrix_allcycles1_mi$time.followup.death.12mnth <- ifelse(matrix_allcycles1_mi$death.dummy.12mnth == 1, matrix_allcycles1_mi$mindeath.index - matrix_allcycles1_mi$mi.index, 12)
total_deaths_mi_12mnth <- sum(matrix_allcycles1_mi$death.dummy.12mnth)
totalyears_deaths_mi_12mnth <- sum(matrix_allcycles1_mi$time.followup.death.12mnth) / 12

mideath_12mnth_rate <- (total_deaths_mi_12mnth / totalyears_deaths_mi_12mnth) * 100

# 5 years (60 months)
matrix_allcycles1_mi$death.dummy.60mnth <- ifelse(matrix_allcycles1_mi$mindeath.index <= (matrix_allcycles1_mi$mi.index + 59), 1, 0)
matrix_allcycles1_mi$time.followup.death.60mnth <- ifelse(matrix_allcycles1_mi$death.dummy.60mnth == 1, matrix_allcycles1_mi$mindeath.index - matrix_allcycles1_mi$mi.index, 60)
total_deaths_mi_60mnth <- sum(matrix_allcycles1_mi$death.dummy.60mnth)
totalyears_deaths_mi_60mnth <- sum(matrix_allcycles1_mi$time.followup.death.60mnth) / 12

mideath_60mnth_rate <- (total_deaths_mi_60mnth / totalyears_deaths_mi_60mnth) * 100 * 5

# 8 years (96 months)
matrix_allcycles1_mi$death.dummy.96mnth <- ifelse(matrix_allcycles1_mi$mindeath.index <= (matrix_allcycles1_mi$mi.index + 95), 1, 0)
matrix_allcycles1_mi$time.followup.death.96mnth <- ifelse(matrix_allcycles1_mi$death.dummy.96mnth == 1, matrix_allcycles1_mi$mindeath.index - matrix_allcycles1_mi$mi.index, 96)
total_deaths_mi_96mnth <- sum(matrix_allcycles1_mi$death.dummy.96mnth)
totalyears_deaths_mi_96mnth <- sum(matrix_allcycles1_mi$time.followup.death.96mnth) / 12

mideath_96mnth_rate <- (total_deaths_mi_96mnth / totalyears_deaths_mi_96mnth) * 100 * 8

# Subsetting those individuals that had MI event and survived one year
matrix_allcycles1_mi12 <- matrix_allcycles1[!is.infinite(matrix_allcycles1$mi12.index), ]

matrix_allcycles1_mi12$mindeath.index <- pmin(matrix_allcycles1_mi12$mideath.index, matrix_allcycles1_mi12$strokedeath.index, matrix_allcycles1_mi12$hfdeath.index,
                                              matrix_allcycles1_mi12$othercvddeath.index, matrix_allcycles1_mi12$chroniccvddeath.index, matrix_allcycles1_mi12$noncvdeath.index,
                                              matrix_allcycles1_mi12$death.index)

# Looking at time to death in 8 years (96 months)
matrix_allcycles1_mi12$death.dummy.96mnth <- ifelse(matrix_allcycles1_mi12$mindeath.index <= (matrix_allcycles1_mi12$mi12.index + 95), 1, 0)

matrix_allcycles1_mi12$time.followup.death.96mnth <- ifelse(matrix_allcycles1_mi12$death.dummy.96mnth == 1, matrix_allcycles1_mi12$mindeath.index - matrix_allcycles1_mi12$mi12.index, 96)

total_deaths_mi12_96mnth <- sum(matrix_allcycles1_mi12$death.dummy.96mnth)
totalyears_deaths_mi12_96mnth <- sum(matrix_allcycles1_mi12$time.followup.death.96mnth) / 12

mideath12_96mnth_rate <- (total_deaths_mi12_96mnth / totalyears_deaths_mi12_96mnth) * 100 * 8

# Looking at time to death in those who had a stroke and survived at least one month
# Subsetting to the requisite dataframe
matrix_allcycles1_stroke2 <- matrix_allcycles1[!is.infinite(matrix_allcycles1$stroke2.index), ]

# Creating a minimum death column
matrix_allcycles1_stroke2$mindeath.index <- pmin(matrix_allcycles1_stroke2$mideath.index, matrix_allcycles1_stroke2$strokedeath.index, matrix_allcycles1_stroke2$hfdeath.index,
                                                 matrix_allcycles1_stroke2$othercvddeath.index, matrix_allcycles1_stroke2$chroniccvddeath.index, matrix_allcycles1_stroke2$noncvdeath.index,
                                                 matrix_allcycles1_stroke2$death.index)
# 12 years (144 months)
matrix_allcycles1_stroke2$death.dummy.144mnth <- ifelse(matrix_allcycles1_stroke2$mindeath.index <= (matrix_allcycles1_stroke2$stroke2.index + 143), 1, 0)
matrix_allcycles1_stroke2$time.followup.death.144mnth <- ifelse(matrix_allcycles1_stroke2$death.dummy.144mnth == 1, matrix_allcycles1_stroke2$mindeath.index - matrix_allcycles1_stroke2$stroke2.index, 144)
total_deaths_stroke2_144mnth <- sum(matrix_allcycles1_stroke2$death.dummy.144mnth)
totalyears_deaths_stroke2_144mnth <- sum(matrix_allcycles1_stroke2$time.followup.death.144mnth) / 12

stroke2death_144mnth_rate <- (total_deaths_stroke2_144mnth / totalyears_deaths_stroke2_144mnth) * 100

# Combining all the rates
df_externalvalid <- data.frame(mnthly_rates = c(
  mideath_12mnth_rate, mideath_60mnth_rate, mideath_96mnth_rate,
  mideath12_96mnth_rate, stroke2death_144mnth_rate
))

rates_names <- c(
  "Mortality rates per 100 person-years after MI (1 year)",
  "Mortality rates per 500 person-years after MI (5 years)",
  "Mortality rates per 800 person-years after MI (8 years)",
  "Mortality rates per 800 person-years after one year post-MI (8 years)",
  "Mortality rates per 100 person-years over 12-year period after 30 days post-Stroke"
)

rates_type_s <- data.frame(type = c(rep("Simulated", 5)))
df_externalvalid1 <- cbind(df_externalvalid, rates_names, rates_type_s)

# Adding the validation targets
df_validtargets <- data.frame(mnthly_rates = c(24.00, 51.00, 65.00, 59.10, 8.23))
rates_type_o <- data.frame(type = c(rep("Observed", 5)))

df_externalvalid2 <- cbind(df_validtargets, rates_names, rates_type_o)

# Combining the simulated and observed dataframes
df_externalvalid3 <- rbind(df_externalvalid1, df_externalvalid2)

df_externalvalid3$rates_names <- factor(df_externalvalid3$rates_names,
  levels = c(
    "Mortality rates per 100 person-years after MI (1 year)",
    "Mortality rates per 500 person-years after MI (5 years)",
    "Mortality rates per 800 person-years after MI (8 years)",
    "Mortality rates per 800 person-years after one year post-MI (8 years)",
    "Mortality rates per 100 person-years over 12-year period after 30 days post-Stroke"
  )
)

# Cleveland dot plot
valid_plot <- ggplot(data = df_externalvalid3, aes(x = mnthly_rates, y = rates_names, color = type, label = format(round(mnthly_rates, digits = 1)))) +
  geom_point(size = 3) +
  geom_line(aes(group = rates_names)) +
  theme_classic() +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  geom_text_repel(point.size = NA, vjust = -1, show.legend = FALSE) +
  xlab("Rates") +
  ylab("Validation targets") +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )

ggsave(
  filename = here("Outputs", paste0("externalvalid_plots_", Sys.Date(), ".png")),
  plot = valid_plot, width = 12, height = 8
)

ggsave(
  filename = here("Outputs", paste0("externalvalid_plots_", Sys.Date(), ".svg")),
  plot = valid_plot, width = 12, height = 10
)
