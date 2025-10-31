library(tidyverse)
library(sybil)
library(BacArena)
library(data.table)
library(wrapr)

validation_sims <- as_named_list(sim_OMM12_aam, sim_OMM12_m10,sim_OMM12_mouse)


##extract SCFA data from models for t4
SCFAs <- c("cpd00029","cpd00141", "cpd00211", "cpd00159", "cpd00221")
SCFA_names <- c("Acetate", "Propionate", "Butyrate", "L-Lactate", "D-Lactate")
SCFA_t4 <- data.frame(matrix(ncol=5, nrow=length(validation_sims)))
colnames(SCFA_t4) <- SCFA_names
sims <- c(1:3)
for (j in 1:length(validation_sims)) {
  sim <- validation_sims[[j]]
  sims[j] <- names(validation_sims)[j]
  
  for (i in 1:length(SCFAs)) {
    ex <- paste("EX", SCFAs[i],"e0", sep = "_")
    value <- sim@medlist[[4]][[ex]][[1]]
    if (is.null(value)==FALSE)
      SCFA_t4[j,i] <- value
    else 
      SCFA_t4[j,i] <- "0"
  }}
rownames(SCFA_t4) <- sims
rm(sim, ex, value)

#add up lactate signals since metabolomic data does not discriminate between D- and L-lactate
SCFA_t4 <- mutate(SCFA_t4,Lactate = `L-Lactate` + `D-Lactate`)

##extract SCFA data from models for t1
SCFAs <- c("cpd00029","cpd00141", "cpd00211", "cpd00159", "cpd00221")
SCFA_names <- c("Acetate", "Propionate", "Butyrate", "L-Lactate", "D-Lactate")
SCFA_t1 <- data.frame(matrix(ncol=5, nrow=length(validation_sims)))
colnames(SCFA_t1) <- SCFA_names
sims <- c(1:3)
for (j in 1:length(validation_sims)) {
  sim <- validation_sims[[j]]
  sims[j] <- names(validation_sims)[j]
  
  for (i in 1:length(SCFAs)) {
    ex <- paste("EX", SCFAs[i],"e0", sep = "_")
    value <- sim@medlist[[1]][[ex]][[1]]
    if (is.null(value)==FALSE)
      SCFA_t1[j,i] <- value
    else 
      SCFA_t1[j,i] <- "0"
  }}
rownames(SCFA_t1) <- sims
rm(sim, ex, value)

#add up lactate signals since metabolomic data does not discriminate between D- and L-lactate
SCFA_t1 <- mutate(SCFA_t1,Lactate = `L-Lactate` + `D-Lactate`)


##extract AA data from models for t4
AAs <- c("cpd00054","cpd00051", "cpd00053", "cpd00023", "cpd00119","cpd00322", "cpd00107", "cpd00039", "cpd00060", "cpd00066", "cpd00129", "cpd00161", "cpd00065", "cpd00156", "cpd00035", "cpd00069", "cpd00041", "cpd00132")
AA_names <- c("Serine",	"Arginine",	"Glutamine", "Glutamic acid", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",	"Phenylalanine", "Proline",	"Threonine", "Tryptophane",	"Valine",	"Alanine", "Tyrosine", "Aspartic acid", "Asparagine")
AA_t4 <- data.frame(matrix(ncol=18, nrow=length(validation_sims)))
colnames(AA_t4) <- AA_names
sims <- c(1:3)
for (j in 1:length(validation_sims)) {
  sim <- validation_sims[[j]]
  sims[j] <- names(validation_sims)[j]
  
  for (i in 1:length(AAs)) {
    ex <- paste("EX", AAs[i],"e0", sep = "_")
    value <- sim@medlist[[4]][[ex]][[1]]
    if (is.null(value)==FALSE)
      AA_t4[j,i] <- value
    else 
      AA_t4[j,i] <- "0"
  }}
rownames(AA_t4) <- sims
rm(sim, ex, value)

#add up leu and ile signals since metabolomic data does not discriminate between leucine and isoleucine
AA_t4 <- mutate(AA_t4,`Isoleucine/Leucine` = `Leucine` + `Isoleucine`)

##extract AA data from models for t1
AAs <- c("cpd00054","cpd00051", "cpd00053", "cpd00023", "cpd00119","cpd00322", "cpd00107", "cpd00039", "cpd00060", "cpd00066", "cpd00129", "cpd00161", "cpd00065", "cpd00156", "cpd00035", "cpd00069", "cpd00041", "cpd00132")
AA_names <- c("Serine",	"Arginine",	"Glutamine", "Glutamic acid", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",	"Phenylalanine", "Proline",	"Threonine", "Tryptophane",	"Valine",	"Alanine", "Tyrosine", "Aspartic acid", "Asparagine")
AA_t1 <- data.frame(matrix(ncol=18, nrow=length(validation_sims)))
colnames(AA_t1) <- AA_names
sims <- c(1:3)
for (j in 1:length(validation_sims)) {
  sim <- validation_sims[[j]]
  sims[j] <- names(validation_sims)[j]
  
  for (i in 1:length(AAs)) {
    ex <- paste("EX", AAs[i],"e0", sep = "_")
    value <- sim@medlist[[1]][[ex]][[1]]
    if (is.null(value)==FALSE)
      AA_t1[j,i] <- value
    else 
      AA_t1[j,i] <- "0"
  }}
rownames(AA_t1) <- sims
rm(sim, ex, value)

#add up leu and ile signals since metabolomic data does not discriminate between leucine and isoleucine
AA_t1 <- mutate(AA_t1,`Isoleucine/Leucine` = `Leucine` + `Isoleucine`)

#calculate differences
SCFA_difference <- SCFA_t4 %>%
  mutate(across(everything(), ~ . - SCFA_t1[[cur_column()]]))

AA_difference <- AA_t4 %>%
  mutate(across(everything(), ~ . - AA_t1[[cur_column()]]))




write.csv(SCFA_difference, "SCFA_difference.csv")
write.csv(AA_difference, "AA_difference.csv")
