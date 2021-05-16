# Script to model hospitalizations and deaths and interface with SyncroSim

library(rsyncrosim)
library(dplyr)
library(purrr)
library(lubridate)

# Setup ----
transformerName <- "Death + Hospitalization Model: In Hospital"

## Connect to SyncroSim ----
myScenario <- scenario()

pipeline <- datasheet(myScenario, "core_Pipeline", lookupsAsFactors = F)
inputData <- datasheet(myScenario, name = "epi_DataSummary", lookupsAsFactors = F)
runSettings <- datasheet(myScenario, "epiModelDeathHospital_RunSettingsInHospital", lookupsAsFactors = F, optional = T)

## Parse settings ----
jurisdictions <- inputData %>%
  filter(TransformerID == runSettings$CaseSourceTransformerID) %>%
  pull(Jurisdiction) %>%
  unique

# Use last day of hospital data as first day of projection
# TODO: handle no hospital data
if(is.na(runSettings$MinimumTimestep)){
  runSettings$MinimumTimestep <- inputData %>%
    filter(Variable == "In Hospital - Daily") %>%
    pull(Timestep) %>%
    sort %>%
    tail(1)  
}


# Project 28 days forward from last hospital date if no projection end date is provided
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
  
  # Get the needed hospital data
  jurisdictionInHospital <- inputData %>%
    filter(Jurisdiction == jurisdictions[j], Variable == "In Hospital - Daily", TransformerID == runSettings$InHospitalSourceTransformerID)
  
  for (i in runSettings$MinimumIteration:runSettings$MaximumIteration){
    #i = 1
    shape <- ((runSettings$HospitalizationRate)^2)/((runSettings$HospitalizationRateSD)^2)
    scale <- ((runSettings$HospitalizationRateSD)^2) / runSettings$HospitalizationRate
    HospitalizationRate <- rgamma(1, shape = shape, scale = scale)
    iterationCases <- jurisdictionCases %>% filter(Iteration == i)
    startInHospital <- jurisdictionInHospital %>% 
      filter(Iteration == i, Timestep == runSettings$MinimumTimestep)
    if(nrow(startInHospital)==0){
      startInHospital <- jurisdictionInHospital %>% 
        filter(is.na(Iteration), Timestep == runSettings$MinimumTimestep)
    }
    if(nrow(startInHospital)==0){
      startInHospital <-0
    } else {
      startInHospital <- startInHospital %>% 
        pull(Value) %>% 
        sort %>%
        tail(1)
    }
    
    myOutput <- addRow(myOutput, list(transformerName[1], 
                          i, as.character(runSettings$MinimumTimestep), 
                          "In Hospital - Daily", jurisdictions[j],
                          NA, NA, NA, startInHospital))
    inHospital <- startInHospital
    for(timestep in (as.integer(as.Date(runSettings$MinimumTimestep))+1):as.integer(as.Date(runSettings$MaximumTimestep))){
      #timestep <- (as.integer(as.Date(runSettings$MinimumTimestep))+1)
      caseNumber <- jurisdictionCases %>% 
        filter(Timestep == as_date(timestep - runSettings$HospitalizationDelay), 
               Iteration == i) %>% 
        pull(Value) %>% sort %>%
        tail(1)
      newInHospital <- rbinom(1,round(caseNumber, digits = 0),HospitalizationRate)
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
        releases <- rbinom(1,round(caseReleases, digits = 0),HospitalizationRate)
      }
      
      inHospital <- inHospital + newInHospital - releases
      if(inHospital < 0)
        inHospital <- 0
      
      myOutput <- addRow(myOutput, list(transformerName[1], 
                                        i, as.character(as_date(timestep)), 
                                        "In Hospital - Daily", jurisdictions[j],
                                        NA, NA, NA, inHospital))
    }
  }
}


# Save run settings back to SyncroSim
saveDatasheet(myScenario, runSettings, "epiModelDeathHospital_RunSettingsInHospital")
saveDatasheet(myScenario, myOutput, "epi_DataSummary", append = T)
