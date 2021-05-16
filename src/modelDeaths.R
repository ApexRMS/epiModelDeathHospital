# Script to model hospitalizations and deaths and interface with SyncroSim

library(rsyncrosim)
library(dplyr)
library(purrr)
library(lubridate)

# Setup ----
transformerName <- "Death + Hospitalization Model: Deaths"

## Connect to SyncroSim ----
myScenario <- scenario()

pipeline <- datasheet(myScenario, "core_Pipeline", lookupsAsFactors = F)
inputData <- datasheet(myScenario, name = "epi_DataSummary", lookupsAsFactors = F)
runSettings <- datasheet(myScenario, "epiModelDeathHospital_RunSettingsDeaths", lookupsAsFactors = F, optional = T)

## Parse settings ----
jurisdictions <- inputData %>%
  filter(TransformerID == runSettings$CaseSourceTransformerID) %>%
  pull(Jurisdiction) %>%
  unique

# Use last day of death data as first day of projection
# TODO: handle no death data
if(is.na(runSettings$MinimumTimestep)){
  runSettings$MinimumTimestep <- inputData %>%
    filter(Variable == "Deaths - Cumulative",
           TransformerID == runSettings$DeathSourceTransformerID) %>%
    pull(Timestep) %>%
    sort %>%
    tail(1)  
}


# Project 28 days forward from last death date if no projection end date is provided
if(is.na(runSettings$MaximumTimestep))
  runSettings$MaximumTimestep <-  runSettings$MinimumTimestep %>%
    ymd %>%
    `+`(28) %>%
    format("%Y-%m-%d")

myOutput <- datasheet(myScenario, "epi_DataSummary", 
                      empty = T, optional = T, 
                      lookupsAsFactors = F)
myOutput$Timestep = as.character(myOutput$Timestep)

envBeginSimulation(length(jurisdictions) * 
                     (runSettings$MaximumIteration - runSettings$MinimumIteration +1) * 
                     (as.integer(as.Date(runSettings$MaximumTimestep)) - 
                        as.integer(as.Date(runSettings$MinimumTimestep)) +1))

for (j in 1:length(jurisdictions)) {

  #j = 1
  # Get the needed case data
  jurisdictionCases <- inputData %>%
    filter(Jurisdiction == jurisdictions[j], Variable == "Cases - Daily", TransformerID == runSettings$CaseSourceTransformerID)
  
  # Get the needed death data
  jurisdictionCumulativeDeaths <- inputData %>%
    filter(Jurisdiction == jurisdictions[j], Variable == "Deaths - Cumulative", TransformerID == runSettings$DeathSourceTransformerID)
  jurisdictionDailyDeaths <- inputData %>%
    filter(Jurisdiction == jurisdictions[j], Variable == "Deaths - Daily", TransformerID == runSettings$DeathSourceTransformerID)
  
  for (i in runSettings$MinimumIteration:runSettings$MaximumIteration){
    #i = 1
    shape <- ((runSettings$CaseFatalityRate)^2)/((runSettings$CaseFatalityRateSD)^2)
    scale <- ((runSettings$CaseFatalityRateSD)^2) / runSettings$CaseFatalityRate
    fatalityRate <- rgamma(1, shape = shape, scale = scale)
    iterationCases <- jurisdictionCases %>% filter(Iteration == i)
    startCumulativeDeaths <- jurisdictionCumulativeDeaths %>% 
      filter(Iteration == i, Timestep == runSettings$MinimumTimestep)
    if(nrow(startCumulativeDeaths)==0){
      startCumulativeDeaths <- jurisdictionCumulativeDeaths %>% 
        filter(is.na(Iteration), Timestep == runSettings$MinimumTimestep)
    }
    if(nrow(startCumulativeDeaths)==0){
      startCumulativeDeaths <-0
    } else {
      startCumulativeDeaths <- startCumulativeDeaths %>% 
        pull(Value) %>% 
        sort %>%
        tail(1)
    }
    
    startDailyDeaths <- jurisdictionDailyDeaths %>% 
      filter(Iteration == i, Timestep == runSettings$MinimumTimestep)
    if(nrow(startDailyDeaths)==0){
      startDailyDeaths <- jurisdictionDailyDeaths %>% 
        filter(is.na(Iteration), Timestep == runSettings$MinimumTimestep)
    }
    if(nrow(startDailyDeaths)==0){
      startDailyDeaths <-0
    } else {
      startDailyDeaths <- startDailyDeaths %>% 
        pull(Value) %>% 
        sort %>%
        tail(1)
    }
    myOutput <- addRow(myOutput, list(transformerName[1], 
                          i, as.character(runSettings$MinimumTimestep), 
                          "Deaths - Cumulative", jurisdictions[j],
                          NA, NA, NA, startCumulativeDeaths))
    if(startDailyDeaths>0)
      myOutput <- addRow(myOutput, list(transformerName[1], 
                                        i, as.character(runSettings$MinimumTimestep), 
                                        "Deaths - Daily", jurisdictions[j],
                                        NA, NA, NA, startDailyDeaths))
    cumulativeDeaths <- startCumulativeDeaths
    
    for(timestep in (as.integer(as.Date(runSettings$MinimumTimestep))+1):as.integer(as.Date(runSettings$MaximumTimestep))){
      #timestep <- (as.integer(as.Date(runSettings$MinimumTimestep))+1)
      caseNumber <- jurisdictionCases %>% 
        filter(Timestep == as_date(timestep - runSettings$DeathDelay), 
               Iteration == i) %>% 
        pull(Value) %>% sort %>%
        tail(1)
      deathNumber <- rbinom(1,round(caseNumber, digits = 0),fatalityRate)
      cumulativeDeaths <- cumulativeDeaths + deathNumber
      myOutput <- addRow(myOutput, list(transformerName[1], 
                                        i, as.character(as_date(timestep)), 
                                        "Deaths - Daily", jurisdictions[j],
                                        NA, NA, NA, deathNumber))
      myOutput <- addRow(myOutput, list(transformerName[1], 
                                        i, as.character(as_date(timestep)), 
                                                        "Deaths - Cumulative", jurisdictions[j],
                                                        NA, NA, NA, cumulativeDeaths))
      envStepSimulation()
    }
  }
}

envEndSimulation()


# Save run settings back to SyncroSim
saveDatasheet(myScenario, runSettings, "epiModelDeathHospital_RunSettingsDeaths")
saveDatasheet(myScenario, myOutput, "epi_DataSummary", append = T)
