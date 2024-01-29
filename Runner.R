## ----Multi-level Regression and Poststratification, 2019 Gen Election----
# Benjamin Mellor
# 21/01/24

## ----Setup, include=FALSE, results="hide", warning=FALSE-----------------
getwd()
setwd('/Volumes/BENSTON/github_repos/MRP_2019Gen')
rm(list=ls(all=TRUE))

## ----installlibs----------------------------------------------------------
#install.packages(
#   "arm", "foreign", "lme4", "haven", "survey", "stanarm", "ggplot2", 
#   "bayesplot", "parallel", "mlogit", "brms"
#       )

## ----loadlibs-------------------------------------------------------------
library("haven")
library("arm")
library("foreign")
library("lme4")
library("survey")
library("dplyr")
library("rstanarm")
library("ggplot2")
library("bayesplot")
library("parallel")
library("tidyverse")
library("mlogit")
library("brms")
theme_set(bayesplot::theme_default())



## ----loadONS--------------------------------------------------------------
# CLEANCLEANCLEAN

age_levels <- c("1", "2", "3", "4", "5", "6")

# read in 
census_data <- read.csv("Data/ONS.csv")
census_data <- census_data %>%
  dplyr::select(
    GSSCode, SexCode, AgeCode, EduCode, Observation
                ) %>%
  rename(
    Sex = SexCode,
    EdLevel = EduCode
  ) %>%
  mutate(AgeCode = as.factor(AgeCode), 
         EdLevel = as.factor(EdLevel),
         Sex = as.factor(Sex)
         ) %>%
  filter(AgeCode != 1, EdLevel != -8)
  
    

# weight
ukpop <- sum(census_data$Observation)
census_data$weight <- census_data$Observation / ukpop

ps <- census_data

## ----loadBESPre2019--------------------------------------------------------
# CLEANCLEANCLEAN

Sample <- read_dta("Data/bes_w18_2019.dta")
Sample <- Sample %>% 
  dplyr::select(
    generalElectionVote, gender, age, ageGroup, pcon_code, p_edlevel
    ) %>% 
  rename(
    Vote = generalElectionVote,
    GSSCode = pcon_code,
    Age = age,
    Sex = gender,
    EdLevel = p_edlevel
  ) %>%
  filter(Vote != 9999)

# Set up for ONS Congruency
age_breaks <- c(-Inf, 15, 24, 34, 49, 64, Inf)
age_labels_func <- c("1", "2", "3", "4", "5", "6")
age_labels <- c("0-15",
                "16-24",
                "25-34",
                "35-49",
                "50-64",
                "65+")

Sample$AgeQ <- as.numeric(
  as.character(Sample$Age)
)
Sample$AgeCode <- cut(
  Sample$AgeQ,
  breaks = age_breaks,
  labels = age_labels_func
)

Sample <- Sample %>%
  select(-ageGroup, -AgeQ, -Age) %>%
  mutate(AgeCode = as.factor(AgeCode), 
         EdLevel = as.factor(EdLevel),
         Sex = as.factor(Sex)
  )

voter_sample <- Sample[complete.cases(Sample),]

## ----loadGEResults---------------------------------------------------------
# CLEANCLEANCLEAN

GE_Results <- read_dta("Data/bes_results_2019.dta")
GE_Results <- GE_Results %>%
  dplyr::select(
    ONSConstID, Winner19, Con19, Lab19, LD19, 
    SNP19, PC19, UKIP19, Green19, Brexit19,
    Other19
  ) %>% 
  rename(
    Winner = Winner19,
    GSSCode = ONSConstID,
  )

## ----loadBESRPS------------------------------------------------------------
# To calculate demographic turnout rates

# read in
BES_rps <- read_dta("Data/bes_rps_2019.dta")
  

BES_rps <- BES_rps %>% 
  dplyr::select(Constit_Code, b01, Q24_CSES, y09, edlevel) %>%
  rename(
    GSSCode = Constit_Code,
    VoterTurnout = b01,
    Age = Q24_CSES,
    Sex = y09,
    EdLevel = edlevel
  ) %>%
  dplyr::filter(Sex %in% c(1, 2)) %>%
  dplyr::filter(VoterTurnout %in% c(1, 2)) %>%
  dplyr::mutate(VoterTurnout = as.factor(VoterTurnout))
  

# categorise and clean for agreement with ONS
BES_rps$AgeQ <- as.numeric(
                  as.character(BES_rps$Age)
                  )
BES_rps$AgeCode <- cut(
                  BES_rps$AgeQ,
                  breaks = age_breaks,
                  labels = age_labels_func
                  )
BES_rps <- BES_rps %>%
  select(-AgeQ, -Age) %>%
  dplyr::mutate(AgeCode = as.factor(AgeCode))

turnout_sample <- BES_rps

## ----SetCores-------------------------------------------------------------
# CLEANCLEANCLEAN

# For parallel QMC
options(mc.cores = parallel::detectCores() - 1)
my_chains <- ifelse(parallel::detectCores() == 1,
                    2,
                    parallel::detectCores() - 1)


## ----MRP-----------------------------------------------------------------

## ----1. EstTurnoutRates--------------------------------------------------

# Fit turnout to sample
turnout_fit <- stan_glmer(
  VoterTurnout ~
      AgeCode + EdLevel + Sex + (1|GSSCode),
  data = turnout_sample,
  family = stats::binomial()
)

# Calculate Constituency Turnout P
const_turnout <- posterior_epred(turnout_fit,
                                 newdata = ps)

const_turnout_P <- colMeans(const_turnout)

# ~ 3GB large
rm(const_turnout)

## ----2. EstVoteOutcome--------------------------------------------------

# Multi Level Regression
voter_fit_polyn <- brm(
  Vote ~
    AgeCode + EdLevel + Sex + (1|GSSCode),
  data = voter_sample,
  family = "categorical",
  seed = 20220302
    )

mse <- brms::predictive_error(fit1, draws = 100) ^ 2

iters <- voter_fit$stanfit@stan_args[[1]]$iter

draws_ps <- brms::posterior_epred(voter_fit,
                                     draws = ifelse(iters > 500,
                                                    500,
                                                    iters),
                                     newdata = ps)

## ----CalcOutcomes--------------------------------------------------

ps$TurnoutP <- const_turnout_P
ps$WeightedTurnout <- ps$TurnoutP * ps$weight

GE_predictions <- draws_ps %>%
  group_by(.draw, .category, GSSCode) %>%
  mutate(share = exp(Vote) / sum(exp(Vote))) %>%
  group_by(AgeCode, Sex, EdLevel, .draw) %>%
  summarize(share = mean(share))


