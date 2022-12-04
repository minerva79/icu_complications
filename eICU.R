#library(tidyverse)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
eICU_patient_data <- read.csv("D:/Downloads/root/eicu-code/data/physionet.org/files/eicu-crd/2.0/patient.csv.gz")
#number of distinct pts, n = 139367
eICU_patient_data %>% count(uniquepid)

head(eICU_patient_data)

#isolate diagnosis codes

eICU_diagnosis_data <- read.csv("D:/Downloads/root/eicu-code/data/physionet.org/files/eicu-crd/2.0/diagnosis.csv.gz")

head(eICU_diagnosis_data)

#admission diagnosis
admission_dx_data <- read.csv("D:/Downloads/root/eicu-code/data/physionet.org/files/eicu-crd/2.0/admissionDx.csv.gz")
ami <- admission_dx_data %>% filter(
  grepl("myocardial", admitdxname)
)

eICU_patient_data_los <-eICU_patient_data %>% 
  select("patientunitstayid","hospitaladmittime24","hospitaladmitoffset","hospitaldischargetime24","hospitaldischargeoffset","uniquepid")

#los of whole admission episode
ami_los <- ami %>% left_join(
  eICU_patient_data_los 
) %>% 
  mutate(
    "LoS (days)" = (`hospitaldischargeoffset` - `hospitaladmitoffset`)/60/24
  )

#n=6933, bassed on admission diag
ami_los %>% count(`uniquepid`)
#removing illogical result
ami_los2 <- ami_los %>% filter(!`LoS (days)` <0)

#finding number of unique patients, n = 6931
ami_los3 <- ami_los2 %>% distinct(`uniquepid`)

#filtering to more than equal 5 days
ami_los_5days <- ami_los2 %>% filter(`LoS (days)` >= 5)

#finding out how many distinct patients with los >= 5, n = 1835
ami_los_5days2 <- ami_los2 %>% filter(`LoS (days)` >= 5) %>% count(`uniquepid`) %>% filter(n == 1)
#checking to make sure number of distinct patients is correct
ami_los_5days3 <- ami_los_5days %>% filter(`uniquepid` %in% ami_los_5days2$uniquepid)


admitdxdate <- ami_los_5days3 %>% select("patientunitstayid","admitdxenteredoffset")

#lab test
lab_data <- read.csv("D:/Downloads/root/eicu-code/data/physionet.org/files/eicu-crd/2.0/lab.csv.gz")
lab_data2 <- lab_data %>% filter(grepl("troponin|calcium",labname,ignore.case = TRUE))
lab_data3 <- lab_data2 %>% filter(`patientunitstayid` %in% ami_los_5days3$patientunitstayid)
lab_data4 <-  lab_data3 %>% arrange(labresultoffset) %>% 
  group_by(`patientunitstayid`,`labname`) %>% mutate("sequence" = row_number()) %>% mutate(labresultoffset = labresultoffset/60/24)

#putting unitstayid back into diagnosis to find possible complication
complication <- eICU_diagnosis_data %>% filter(`patientunitstayid` %in% ami_los_5days3$patientunitstayid) %>% left_join(
  admitdxdate
) %>% filter(icd9code != "") %>% 
  mutate(
    "renal" = case_when(
      grepl("584.9|N17.9", icd9code) ~ 1,
      T ~ 0
    )
  ) %>% filter(`renal` == 1)

#lab data for renal fn
lab_data_renal <- lab_data %>% filter(grepl("creatinine|ammonia|gfr|bun",labname,ignore.case = TRUE))
lab_data_renal2 <- lab_data_renal %>% filter(`patientunitstayid` %in% ami_los_5days3$patientunitstayid) %>% 
  # arrange(labresultoffset) %>% 
  group_by(`patientunitstayid`,`labresultoffset`,`labname`) %>% summarise(
    labresult = mean(labresult)
  ) %>% mutate(`complication`= case_when(
    `patientunitstayid` %in% complication$patientunitstayid ~ 1,
    T ~ 0
  ))
df <- lab_data_renal2 %>% ungroup() %>% group_by(`labresultoffset`,`labname`,`complication`) %>% 
  summarise(
    labresult = mean(labresult)) %>% filter(labname %in% c("creatinine","BUN"))

#reading data for patients with lung problem
resp_prob <- read.csv("C:/Users/ch000/Downloads/data.csv") %>% filter(`lung_problem` == "True")
resp_no_prob <- read.csv("C:/Users/ch000/Downloads/data.csv") %>% filter(`lung_problem` == "False")

#filtering lab data for respiratory problem
lab_data_resp <- lab_data %>% filter(grepl("HCO3|paco2|pao2",labname,ignore.case = TRUE)) %>% filter(`patientunitstayid` %in% ami_los_5days3$patientunitstayid)
lab_data_resp2 <- lab_data_resp %>% 
  group_by(`patientunitstayid`,`labresultoffset`,`labname`) %>% summarise(
    labresult = mean(labresult)
  ) %>% mutate(`complication`= case_when(
    `patientunitstayid` %in% resp_prob$patientunitstayid ~ as.character("lung_problem"),
    `patientunitstayid` %in% resp_no_prob$patientunitstayid ~ as.character("no_lung_problem")
  ))
df2 <- lab_data_resp2 %>% ungroup() %>% group_by(`labresultoffset`,`labname`,`complication`) %>% 
  summarise(
    labresult = mean(labresult)) %>% filter(grepl("HCO3|paco2|pao2",labname,ignore.case = TRUE))
df2 <- df2 %>%  spread(labname, labresult)
write.csv(lab_data_resp, "C:/Users/ch000/Downloads/HCO3paC02PaO2.csv",row.names = FALSE)

# categorise to standards for each paramater
df2$HCO3_group <- cut(df2$HCO3 , breaks=c(-Inf, 22, 26, Inf))
levels(df2$HCO3_group) <- c("below", "normal", "above")
df2$paCO2_group <- cut(df2$paCO2, breaks=c(-Inf, 35, 45, Inf))
levels(df2$paCO2_group) <- c("below", "normal", "above")
df2$paO2_group <- cut(df2$paO2, breaks=c(-Inf, 75, 100, Inf))
levels(df2$paO2_group) <- c("below", "normal", "above")

# plotting data
# df2$complication <- factor((df2$complication + 0), labels=c("lung_problem", "no_lung_problem")) 
df2 <- df2 %>% mutate(time = labresultoffset/1440)

# ggplot and putting trendline with generalised additive modelling
ggplot(df2 %>% filter(!is.na(complication)), aes(time, HCO3, col=HCO3_group)) +
  labs(y = "HCO3(mmol/L)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=22), linetype = "dashed") +
  geom_hline(aes(yintercept=26), linetype = "dashed") + 
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, paCO2, col=paCO2_group)) +
  labs(y = "paCO2(mm Hg)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=35), linetype = "dashed") +
  geom_hline(aes(yintercept=45), linetype = "dashed") + 
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, paO2, col=paO2_group)) +
  labs(y = "paO2(mm Hg)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=75), linetype = "dashed") +
  geom_hline(aes(yintercept=100), linetype = "dashed") + 
  theme_bw()

#getting the top 20 icd-10 codes
a <- strsplit(complication$icd9code,split = ", ") %>% unlist()
icd_9_df <- as.data.frame(a) %>% count(a) %>% arrange(-n) %>% filter(grepl("^[[:alpha:]]",a))
icd_9_df <- icd_9_df %>% mutate_at("a",str_replace,"\\.","")
icd_data <- read.csv("C:/Users/ch000/Downloads/icd10.csv") %>% select("A000","Cholera.due.to.Vibrio.cholerae.01..biovar.cholerae") %>% 
  rename("a" = "A000")
icd_9_df <- icd_9_df[1:20,] %>% left_join(icd_data) %>% mutate(
  "group" = case_when(
    `a` %in% "J9600" ~ as.character("respiratory"),
    `a` %in% "N179" ~ as.character("renal"),
    grepl("I21|I50|R570|I480|I469|I2510|I959",`a`) ~ as.character("cardiac"),
    T ~ as.character("others")
  )
)

#plotting barchart for top 20 icd-10 codes
ggplot(icd_9_df, aes(x =reorder(a,-n), y = n,fill= group)) +geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.8)) +
  labs(x= "ICD-10 codes", y = "number of patients with respective complications") +
  scale_fill_manual(values=c("red","grey","green","blue"))

icd_9_desc_df <- icd_9_df %>% select("a","Cholera.due.to.Vibrio.cholerae.01..biovar.cholerae") %>% 
  rename("ICD-10 codes"= "a",
         "ICD-10 description" = "Cholera.due.to.Vibrio.cholerae.01..biovar.cholerae")
write.csv(icd_9_desc_df,"C:/Users/ch000/Downloads/top_20icd10.csv",row.names = FALSE )

# ami <- eICU_diagnosis_data %>% 
#   filter(
#     icd9code == "414.00, I25.10"
#   ) %>% View()
# 
# ami <- eICU_diagnosis_data %>% select(icd9code) %>% View()
# 
# ami_2 <- eICU_diagnosis_data %>% 
#   filter(
#     diagnosisoffset == 72
#   ) %>% View()

#lab data for electrolytes
renal_prob <- read.csv("C:/Users/ch000/Downloads/data.csv") %>% filter(`kidney_problems` == "True")
renal_no_prob <- read.csv("C:/Users/ch000/Downloads/data.csv") %>% filter(`kidney_problems` == "False")
lab_data_elec<- lab_data %>% filter(`labname` %in% c("calcium","phosphate","glucose","potassium","magnesium","sodium","chloride"))
lab_data_elec2 <- lab_data_elec %>% filter(`patientunitstayid` %in% ami_los_5days3$patientunitstayid) %>% 
  group_by(`patientunitstayid`,`labresultoffset`,`labname`) %>% summarise(
    labresult = mean(labresult)
  ) %>% mutate(`complication`= case_when(
    `patientunitstayid` %in% c(resp_prob$patientunitstayid,renal_prob$patientunitstayid) ~ as.character("problem"),
    `patientunitstayid` %in% c(resp_no_prob$patientunitstayid,renal_no_prob$patientunitstayid) ~ as.character("no_problem")
  ))
df2 <- lab_data_elec2 %>% ungroup() %>% group_by(`labresultoffset`,`labname`,`complication`) %>% 
  summarise(
    labresult = mean(labresult)) %>% filter(`labname` %in% c("calcium","phosphate","glucose","potassium","magnesium","sodium","chloride"))
df2 <- df2 %>%  spread(labname, labresult)
write.csv(lab_data_elec, "C:/Users/ch000/Downloads/electrolytes.csv",row.names = FALSE)

# categorise to standards for each paramater
df2$calcium_group <- cut(df2$calcium, breaks=c(-Inf, 8.4, 10.2, Inf))
levels(df2$calcium_group) <- c("below", "normal", "above")
df2$phosphate_group <- cut(df2$phosphate, breaks=c(-Inf, 3.0, 4.5, Inf))
levels(df2$phosphate_group) <- c("below", "normal", "above")
df2$glucose_group <- cut(df2$glucose , breaks=c(-Inf, 70, 110, Inf))
levels(df2$glucose_group) <- c("below", "normal", "above")
df2$potassium_group <- cut(df2$potassium, breaks=c(-Inf, 3.5, 5, Inf))
levels(df2$potassium_group) <- c("below", "normal", "above")
df2$magnesium_group <- cut(df2$magnesium, breaks=c(-Inf, 1.7, 2.2, Inf))
levels(df2$magnesium_group) <- c("below", "normal", "above")
df2$sodium_group <- cut(df2$sodium , breaks=c(-Inf, 136, 146, Inf))
levels(df2$sodium_group) <- c("below", "normal", "above")
df2$chloride_group <- cut(df2$chloride, breaks=c(-Inf, 95, 105, Inf))
levels(df2$chloride_group) <- c("below", "normal", "above")

# plotting data
df2 <- df2 %>% mutate(time = labresultoffset/1440)

# ggplot and putting trendline with generalised additive modelling
ggplot(df2 %>% filter(!is.na(complication)), aes(time, calcium, col=calcium_group)) +
  labs(y = "calcium(mg/dL)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=8.4), linetype = "dashed") +
  geom_hline(aes(yintercept=10.2), linetype = "dashed") +
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, phosphate, col=phosphate_group)) +
  labs(y = "phosphate(mg/dL)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=3.0), linetype = "dashed") +
  geom_hline(aes(yintercept=4.5), linetype = "dashed") + 
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, glucose, col=glucose_group)) +
  labs(y = "glucose(mg/dL)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=70), linetype = "dashed") +
  geom_hline(aes(yintercept=110), linetype = "dashed") + 
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, potassium, col=potassium_group)) +
  labs(y = "potassium(mmol/L)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=3.5), linetype = "dashed") +
  geom_hline(aes(yintercept=5), linetype = "dashed") +
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, magnesium, col=magnesium_group)) +
  labs(y = "magnesium(mg/dL)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=1.7), linetype = "dashed") +
  geom_hline(aes(yintercept=2.2), linetype = "dashed") + 
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, sodium, col=sodium_group)) +
  labs(y = "sodium(mmol/L)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=136), linetype = "dashed") +
  geom_hline(aes(yintercept=146), linetype = "dashed") + 
  theme_bw()

ggplot(df2 %>% filter(!is.na(complication)), aes(time, chloride, col=chloride_group)) +
  labs(y = "chloride(mmol/L)", x = "time(days)") +
  geom_point() +
  stat_smooth(method="gam", col="black") +
  facet_wrap(~complication) +
  geom_hline(aes(yintercept=95), linetype = "dashed") +
  geom_hline(aes(yintercept=105), linetype = "dashed") +
  theme_bw()