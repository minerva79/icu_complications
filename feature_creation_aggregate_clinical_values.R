###########################################################################
#### Feature Creation - Aggregated scores for subject's clinical tests ####
#### SHADE '22 Datathon Team 6 / 2 - 4 December 2022                   ####
###########################################################################


# Packages Dependencies ---------------------------------------------------

lib <- c("plyr", "tidyverse", "mgcv")
invisible(lapply(lib, require, character=TRUE))
rm(lib)


# Data source dependencies ------------------------------------------------

# eICU Research Database - lab
lab_df <- read_csv("lab.csv")

# subjects with AMI (n= 1835)
subjects <- read_csv("patientunitstayid.csv") %>% unlist
# comorbidities/ complication diagnosis for subjects
dx_df <- read_csv('dx_data.csv')


# Custom Function ---------------------------------------------------------

compareFunction <- function(i, data, model, parameter="creatinine"){
  ### Function fit subject's lab time to the generalised additive model (model)
  ### and compare the subject's lab value to the fitted value.
  ### - Assigned value of 1, 0 and -1, if lab value is above, within or below 
  ### cohort fitted threshold.
  ### - Return aggregated value (mean); range of -1 to 1.
  
  e101 <- data %>%  
    filter(patientunitstayid == i)
  e102 <- augment(model, newdata = e101) %>% 
    select(.fitted, .se.fit, one_of(parameter))
  names(e102) <- c("fitted", "se", "true_val")
  e102 <- na.exclude(e102)
  e102 <- e102 %>% 
    mutate(score = ifelse(true_val > fitted + se, 1,
                          ifelse(true_val < fitted - se, -1, 0)))
  return(mean(e102$score))
}




# Exploration of temporal renal trends ------------------------------------

# getting the renal panel results from the subjects
subject_renal_panel <- lab_df %>% filter(patientunitstayid %in% subjects)
subject_renal_panel <- subject_renal_panel %>%
  filter(labname %in% c("BUN", "creatinine"))

# subject_renal_panel %>% 
#   group_by(patientunitstayid, labresultoffset) %>% 
#   filter(labname == "BUN") %>%  
#   tally() %>%
#   arrange(desc(n))

# aggregate duplicated subject tests by mean
subject_renal_panel <- subject_renal_panel %>% 
  group_by(patientunitstayid, labresultoffset, labname) %>% 
  summarise_at(vars(labresult), funs(mean))

subject_renal_panel <- subject_renal_panel %>% 
  spread(labname, labresult)

# create BUN to creatinine ratio
subject_renal_panel <- subject_renal_panel %>% 
  mutate(BUN_creatinine_ratio = BUN/creatinine)

# categorise to standards for each paramater
subject_renal_panel$BUN_group <- cut(subject_renal_panel$BUN, breaks=c(-Inf, 7, 20, Inf))
levels(subject_renal_panel$BUN_group) <- c("below", "normal", "above")

subject_renal_panel$creatinine_group <- cut(subject_renal_panel$creatinine, breaks=c(-Inf, .7, 1.2, Inf))
levels(subject_renal_panel$creatinine_group) <- c("below", "normal", "above")

subject_renal_panel$BUN_creatinine_ratio_group <- cut(subject_renal_panel$BUN_creatinine_ratio, breaks=c(-Inf, 12, 20.1, Inf))
levels(subject_renal_panel$BUN_creatinine_ratio_group) <- c("below", "normal", "above")

# join with the diagnosis code
subject_renal_panel <- subject_renal_panel %>% 
  left_join(dx_df)

# plotting data
df <- subject_renal_panel
df$kidney_problems <- factor((df$kidney_problems + 0), labels=c("no kidney dx", "kidney dx")) 
df <- df %>% 
  mutate(time = labresultoffset/1440)

# ggplot and putting trendline with generalised additive modelling
ggplot(df %>% filter(!is.na(kidney_problems) & !is.na(BUN_group)), aes(time, BUN, col=BUN_group)) +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~kidney_problems) +
  ylab("Blood urea nitrogen (BUN) (mg/dL)") +
  xlab("Offset (Days)") +
  guides(col=guide_legend("Clinical Range")) +
  geom_hline(aes(yintercept=7), linetype = "dashed") +
  geom_hline(aes(yintercept=20), linetype = "dashed") + 
  theme_bw()

ggplot(df %>% filter(!is.na(kidney_problems) & !is.na(creatinine)), aes(time, creatinine, col=creatinine_group)) +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~kidney_problems) +
  geom_hline(aes(yintercept=.7), linetype = "dashed") +
  geom_hline(aes(yintercept=1.2), linetype = "dashed") + 
  theme_bw() +
  ylab("Creatinine (mg/dL)") +
  xlab("Offset (Days)") +
  guides(col=guide_legend("Clinical Range"))

ggplot(df %>% filter(!is.na(kidney_problems)), aes(time, BUN_creatinine_ratio, col=BUN_creatinine_ratio_group)) +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~kidney_problems) +
  geom_hline(aes(yintercept=12), linetype = "dashed") +
  geom_hline(aes(yintercept=20.1), linetype = "dashed") + 
  theme_bw() +
  ylab("BUN is to Creatinine ratio") +
  xlab("Offset (Days)") +
  guides(col=guide_legend("Clinical Range"))

# remove the outlier for BUN_creatinine_ratio (point with 600)
ggplot(df %>% filter(!is.na(kidney_problems) & BUN_creatinine_ratio<500), aes(time, BUN_creatinine_ratio, col=BUN_creatinine_ratio_group)) +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~kidney_problems) +
  geom_hline(aes(yintercept=12), linetype = "dashed") +
  geom_hline(aes(yintercept=20.1), linetype = "dashed") + 
  theme_bw() +
  ylab("BUN is to Creatinine ratio") +
  xlab("Offset (Days)") +
  guides(col=guide_legend("Clinical Range"))

# write_csv(df, "subject_renal_panel.csv")





# Comparing parameter to expected value -----------------------------------


# Each parameter model predict the expected mean +/- std error for parameter's value 
# at given time (based on whole icu population).
# Each patient's test result are compared with the modeled value:
#     - if true value > model mean + se, score = 1
#     - if true value < model mean - se, score = -1
#     - if true value is between the mean and se range, score = 0
# All score aggregated to mean.
# Example:
#  Patient 144433, has creatinine of 3.03, 2.55, 2.85 and 3.26 at time (0.426, 1.18, 1.51, 2.19)
#  the modelled values at the same time are 1.6 (+/- different se).
#  Patient creatinine_comparison score is thus 1.

# # A tibble: 4 × 4
# fitted     se true_val score
# <dbl>  <dbl>    <dbl> <dbl>
# 1   1.60 0.0186     3.03     1
# 2   1.60 0.0167     2.55     1
# 3   1.60 0.0168     2.85     1
# 4   1.60 0.0180     3.26     1 


creatinine_mod <- gam(creatinine ~ s(time), data=df)


new_df <- data.frame(patientunitstayid = unique(df$patientunitstayid))
new_df$creatinine_comparison <- sapply(unique(new_df$patientunitstayid), 
            function(x)compareFunction(x, df, creatinine_mod, "creatinine"))

bun_mod <- gam(BUN ~ s(time), data=df)
new_df$bun_comparison <- sapply(unique(new_df$patientunitstayid), 
                                function(x)compareFunction(x, df, bun_mod, "BUN"))

bun_creatinine_mod <- gam(BUN_creatinine_ratio ~ s(time), data=df)
new_df$bun_creatinine_comparison <- sapply(unique(new_df$patientunitstayid), 
                                function(x)compareFunction(x, df, bun_mod, "BUN_creatinine_ratio"))


# write_csv(new_df, "subject_renal_panel_new_feature.csv")


# Modelling the other electrolytes and parameters -------------------------
var <- c("calcium","phosphate","glucose","potassium","magnesium","sodium",
         "chloride", "HCO3", "paCO2", "paO2")
df2 <- lab_df %>% filter(patientunitstayid %in% subjects & labname %in% var)
df2 <- df2 %>% 
  group_by(patientunitstayid, labresultoffset, labname) %>% 
  summarise_at(vars(labresult), funs(mean))
df3 <- df2 %>% spread(labname, labresult)
# saveRDS(df3, "output_data_for_modelling.rds")

df3 <- df3 %>% mutate(time = labresultoffset/1440)

calcium_mod <- gam(calcium ~ s(time), data=df3)
plot(calcium_mod)
new_df2 <- new_df
new_df2$calcium_comparison <- sapply(unique(new_df$patientunitstayid), 
                                           function(x)compareFunction(x, df3, calcium_mod, "calcium"))

hco3_mod <- gam(HCO3 ~ s(time), data=df3)
plot(hco3_mod)
new_df2$hco3_comparison <- sapply(unique(new_df$patientunitstayid), 
                                     function(x)compareFunction(x, df3, hco3_mod, "HCO3"))

# df3 %>% filter(patientunitstayid == 141894) %>% select(HCO3, time)


phosphate_mod <- gam(phosphate ~ s(time), data=df3)
plot(phosphate_mod)
new_df2$phosphate_comparison <- sapply(unique(new_df$patientunitstayid), 
                                  function(x)compareFunction(x, df3, phosphate_mod, "phosphate"))
df3 %>% filter(patientunitstayid == 141894) %>% select(phosphate, time)

# saveRDS(new_df2, "engineered_features.rds")

glucose_mod <- gam(glucose ~ s(time), data=df3)
plot(glucose_mod)
new_df2$glucose_comparison <- sapply(unique(new_df$patientunitstayid), 
                                       function(x)compareFunction(x, df3, glucose_mod, "glucose"))

# saveRDS(new_df2, "engineered_features.rds")
# saveRDS(df3, "working_dat.rds")
# saveRDS(compareFunction, "comparative_f(x).rds")

potassium_mod <- gam(potassium ~ s(time), data=df3)
plot(potassium_mod)
new_df2$potassium_comparison <- sapply(unique(new_df$patientunitstayid), 
                                     function(x)compareFunction(x, df3, potassium_mod, "potassium"))
# saveRDS(new_df2, "engineered_features.rds")
# write_csv(new_df2, "engineered_features.csv")

magnesium_mod <- gam(magnesium ~ s(time), data=df3)
plot(magnesium_mod)
new_df2$magnesium_comparison <- sapply(unique(new_df$patientunitstayid), 
                                       function(x)compareFunction(x, df3, magnesium_mod, "magnesium"))
# saveRDS(new_df2, "engineered_features.rds")
# write_csv(new_df2, "engineered_features.csv")

sodium_mod <- gam(sodium ~ s(time), data=df3)
plot(sodium_mod)
new_df2$sodium_comparison <- sapply(unique(new_df$patientunitstayid), 
                                       function(x)compareFunction(x, df3, sodium_mod, "sodium"))
# saveRDS(new_df2, "engineered_features.rds")
# write_csv(new_df2, "engineered_features.csv")

chloride_mod <- gam(chloride ~ s(time), data=df3)
plot(chloride_mod)
new_df2$chloride_comparison <- sapply(unique(new_df$patientunitstayid), 
                                    function(x)compareFunction(x, df3, chloride_mod, "chloride"))
# saveRDS(new_df2, "engineered_features.rds")
# write_csv(new_df2, "engineered_features.csv")

paCO2_mod <- gam(paCO2 ~ s(time), data=df3)
plot(paCO2_mod)
new_df2$paCO2_comparison <- sapply(unique(new_df$patientunitstayid), 
                                      function(x)compareFunction(x, df3, paCO2_mod, "paCO2"))
# saveRDS(new_df2, "engineered_features.rds")
# write_csv(new_df2, "engineered_features.csv")

paO2_mod <- gam(paO2 ~ s(time), data=df3)
plot(paO2_mod)
new_df2$paO2_comparison <- sapply(unique(new_df$patientunitstayid), 
                                   function(x)compareFunction(x, df3, paO2_mod, "paO2"))
# saveRDS(new_df2, "engineered_features.rds")
write_csv(new_df2, "engineered_features.csv")

# 
# compareFunction(144433, df, creatinine_mod, "creatinine")
# 
# compareFunction(144954, df, creatinine_mod, "creatinine")
# compareFunction(141894, df, creatinine_mod, "creatinine")
# 
# compareFunction(144954, df, bun_mod, "BUN")
# compareFunction(141894, df, bun_mod, "BUN")
# 
# 
# compareFunction(144954, df, bun_creatinine_mod, "BUN_creatinine_ratio")
# compareFunction(141894, df, bun_creatinine_mod, "BUN_creatinine_ratio")


# Initial modeeling codes -----------------------------------------------------
### DO NOT RUN ###

# model_1 <- gam(creatinine ~ s(time), data=df %>% filter(kidney_problems == "kidney dx"))
# model_2 <- gam(creatinine ~ s(time, k=5), data=df %>% filter(kidney_problems == "no kidney dx"))
# plot(model_1)
# plot(model_2)
# 
# plotdat <- augment(model_1, newdata = data.frame(time = seq(-10, 80, length=500)))
# ggplot(plotdat, aes(x = time)) +
#   geom_point(data = df %>% filter(kidney_problems == "kidney dx"), aes(y=creatinine), alpha=.05) +
#   geom_line(aes(y = .fitted), col="black") +
#   geom_line(aes(y = .fitted + .se.fit), linetype="dashed") +
#   geom_line(aes(y = .fitted - .se.fit), linetype="dashed") + 
#   theme_bw()
# 
# plotdat <- augment(model_2, newdata = data.frame(time = seq(-10, 80, length=500)))
# ggplot(plotdat, aes(x = time)) +
#   geom_point(data = df %>% filter(kidney_problems == "no kidney dx"), aes(y=creatinine), alpha=.05) +
#   geom_line(aes(y = .fitted), col="black") +
#   geom_line(aes(y = .fitted + .se.fit), linetype="dashed") +
#   geom_line(aes(y = .fitted - .se.fit), linetype="dashed") + 
#   theme_bw()
# 
# means <- function(x) mean(x, na.rm = TRUE)
# df %>% group_by(kidney_problems) %>% summarise_at(vars(creatinine, BUN, BUN_creatinine_ratio),
#                                                   funs(means))
# 
# m1_fit <- predict(model_1, new_data = plotdat, se.fit=TRUE)$fit
# m1_se <- predict(model_1, new_data = plotdat, se.fit=TRUE)$se.fit