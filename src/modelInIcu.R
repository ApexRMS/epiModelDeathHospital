# Script to model hospitalizations and deaths and interface with SyncroSim

library(rsyncrosim)
library(dplyr)
library(purrr)
library(lubridate)

# Setup ----
transformerName <- "Death + Hospitalization Model: In ICU"

## Connect to SyncroSim ----
myScenario <- scenario()

pipeline <- datasheet(myScenario, "core_Pipeline", lookupsAsFactors = F)
inputData <- datasheet(myScenario, name = "epi_DataSummary", lookupsAsFactors = F)
runSettings <- datasheet(myScenario, "epiModelDeathHospital_RunSettingsInIcu", lookupsAsFactors = F, optional = T)

## Parse settings ----
jurisdictions <- inputData %>%
  filter(TransformerID == runSettings$CaseSourceTransformerID) %>%
  pull(Jurisdiction) %>%
  unique

# Use last day of ICU data as first day of projection
# TODO: handle no ICU data
if(is.na(runSettings$MinimumTimestep)){
  runSettings$MinimumTimestep <- inputData %>%
    filter(Variable == "In ICU - Daily") %>%
    pull(Timestep) %>%
    sort %>%
    tail(1)  
}


# Project 28 days forward from last ICU date if no projection end date is provided
if(is.na(runSettings$MaximumTimestep))
  runSettings$MaximumTimestep <-  runSettings$MinimumTimestep %>%
    ymd %>%
    `+`(28) %>%
    format("%Y-%m-%d")

myOutput <- datasheet(myScenario, "epi_DataSummary", 
                      empty = T, optional = T, 
                      lookupsAsFactors = F)
myOutput$Timestep = as.character(myOutput$Timestep)

for (j in 1:length(jurisdictions)) {

  #j = 1
  # Get the needed case data
  jurisdictionCases <- inputData %>%
    filter(Jurisdiction == jurisdictions[j], Variable == "Cases - Daily", TransformerID == runSettings$CaseSourceTransformerID)
  
  # Get the needed ICU data
  jurisdictionInICU <- inputData %>%
    filter(Jurisdiction == jurisdictions[j], Variable == "In ICU - Daily", TransformerID == runSettings$InICUSourceTransformerID)
  
  for (i in runSettings$MinimumIteration:runSettings$MaximumIteration){
    #i = 1
    shape <- ((runSettings$HospitalizationRate)^2)/((runSettings$HospitalizationRateSD)^2)
    scale <- ((runSettings$HospitalizationRateSD)^2) / runSettings$HospitalizationRate
    HospitalizationRate <- rgamma(1, shape = shape, scale = scale)
    
    shape <- ((runSettings$ProportionICU)^2)/((runSettings$ProportionICUSD)^2)
    scale <- ((runSettings$ProportionICUSD)^2) / runSettings$ProportionICU
    ICUProportion <- rgamma(1, shape = shape, scale = scale)
    
    ICURate <- HospitalizationRate*ICUProportion
    
    iterationCases <- jurisdictionCases %>% filter(Iteration == i)
    startInICU <- jurisdictionInICU %>% 
      filter(Iteration == i, Timestep == runSettings$MinimumTimestep)
    if(nrow(startInICU)==0){
      startInICU <- jurisdictionInICU %>% 
        filter(is.na(Iteration), Timestep == runSettings$MinimumTimestep)
    }
    if(nrow(startInICU)==0){
      startInICU <-0
    } else {
      startInICU <- startInICU %>% 
        pull(Value) %>% 
        sort %>%
        tail(1)
    }
    
    myOutput <- addRow(myOutput, list(transformerName[1], 
                          i, as.character(runSettings$MinimumTimestep), 
                          "In ICU - Daily", jurisdictions[j],
                          NA, NA, NA, startInICU))
    inICU <- startInICU
    for(timestep in (as.integer(as.Date(runSettings$MinimumTimestep))+1):as.integer(as.Date(runSettings$MaximumTimestep))){
      #timestep <- (as.integer(as.Date(runSettings$MinimumTimestep))+1)
      caseNumber <- jurisdictionCases %>% 
        filter(Timestep == as_date(timestep - runSettings$HospitalizationDelay), 
               Iteration == i) %>% 
        pull(Value) %>% sort %>%
        tail(1)
      newInICU <- rbinom(1,round(caseNumber, digits = 0),ICURate)
      caseReleases <- 0
      
      # releases are older hospitalizations
      caseReleases <- jurisdictionCases %>% 
        filter(Timestep == as_date(timestep - runSettings$HospitalizationDelay - runSettings$HospitalizationDelay), 
               Iteration == i) %>% 
        pull(Value) %>% sort %>%
        tail(1)
      
      if(length(caseReleases)==0){
        releases <- 0
      } else {
        releases <- rbinom(1,round(caseReleases, digits = 0),ICURate)
      }
      
      inICU <- inICU + newInICU - releases
      if(inICU < 0){
        inICU <- 0
      }
      
      myOutput <- addRow(myOutput, list(transformerName[1], 
                                        i, as.character(as_date(timestep)), 
                                        "In ICU - Daily", jurisdictions[j],
                                        NA, NA, NA, inICU))
    }
  }
}


# Save run settings back to SyncroSim
saveDatasheet(myScenario, runSettings, "epiModelDeathHospital_RunSettingsInIcu")
saveDatasheet(myScenario, myOutput, "epi_DataSummary", append = T)
