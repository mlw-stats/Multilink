library(lme4)
library(tidyverse)
library(performance)
library(extraDistr)

#A function to simulate data for a cluster randomized trial

simulate_data <- function(n_cluster, cluster_size, nRegions) {
  
  getGammaFromMeanSd<-function(mean,sd){
    a<-mean^2/sd^2
    b<-sd^2/mean
    return(c(a,b))
  }
  
  ## create a hospital level data frame
  if(nRegions==3){
    grH<-expand.grid(paste(sep="","H",1:n_cluster),c("Mw","Tz-1", "Tz-2"))
  }else if(nRegions==2){
    grH<-expand.grid(paste(sep="","H",1:n_cluster),c("Mw","Tz"))
  }
  
  hDat<-data.frame(
    site=paste(sep="_",grH[,2],grH[,1]),
    country=grH[,2],
    studyArm="control",
    propHivScreening=rnorm(n=nrow(grH),mean=0.364,sd=0.025),
    propDiabScreening=rnorm(n=nrow(grH),mean=0.247,sd=0.02),
    propHyperScreening=rnorm(n=nrow(grH),mean=0.119,sd=0.01),
    rateevent=rgamma(n=nrow(grH),shape=getGammaFromMeanSd(mean=0.459,sd=0.1)[1],scale =getGammaFromMeanSd(mean=0.459,sd=0.1)[2])
  )
  
  for(c in unique(hDat$country)){
    idInt<-paste(sep="_",c,paste(sep="","H",sample(1:n_cluster,size=n_cluster/2,replace=F))) # select the intervention hospitals
    hDat$studyArm[hDat$site %in% idInt]<-"intervention"
  }
  
  # simulate the effect of the intervention (as per protocol we assume a drop in readmission or death from 40% to 30.3%)
  hDat$rateevent[hDat$studyArm=="intervention"]<-rgamma(n=sum(hDat$studyArm=="intervention"),shape=getGammaFromMeanSd(mean=0.322,sd=0.1)[1],scale =getGammaFromMeanSd(mean=0.322,sd=0.1)[2])
 
  ## create a patient level data frame
  if(nRegions==3){
    gr<-expand.grid(1:cluster_size,paste(sep="","H",1:n_cluster),c("Mw","Tz-1", "Tz-2"))
  }else if(nRegions==2){
    gr<-expand.grid(1:cluster_size,paste(sep="","H",1:n_cluster),c("Mw","Tz"))
  }
  
  pDat<-data.frame(
    pid=paste(sep="_",gr[,3],gr[,2],gr[,1]), # patient ID
    site=paste(sep="_",gr[,3],gr[,2]) # hospital ID
  )
  
  pDat<-dplyr::left_join(x=pDat,y=hDat,by="site")
  
  pDat<-pDat %>%
    dplyr::mutate(
      date_at_entry=sample(seq(as.Date("2024-02-01"), as.Date("2025-02-28"), by = "days"), n(), replace = TRUE),
      age=18+rpois(n=n(),lambda=10),
      sex=sample(c("M","F"),replace=T,size=n()), # assumed 50-50 proportions
      height= round(ifelse(sex == "F", rnorm(n = sum(sex == "F"), mean=158, sd=5), rnorm(n = sum(sex == "M"), mean = 168, sd=5))),
      weight=round(rnorm(n=n(),mean=20*(height/100)^2,sd=3)),
      hiv=rep(0,nrow(pDat)), #to be derived from elig_hiv1 and samp_hiv
      diabetes=rep(0,nrow(pDat)),
      hypertension=rep(0,nrow(pDat)),
      comorb = NA,
      event_count = rpois(n(), lambda = rateevent),
    ) %>%
    dplyr::mutate(
      ageCentred=age-median(age)
    )
  
  # generate data on comorbidities, guaranteeing that every participant has at least 2 morbidities and deriving the overall comorbidity status
  for(j in 1:nrow(pDat)){
    comorb<-0
    
    while(comorb<2){
      hiv<-rbinom(n=1,size=1,prob=pDat$propHivScreening[j])
      diabetes<-rbinom(n=1,size=1,prob=pDat$propDiabScreening[j])
      hypertension<-rbinom(n=1,size=1,prob=pDat$propHyperScreening[j])
      
      comorb<-hiv+diabetes+hypertension
    }
    
    pDat$hiv[j]<-hiv
    pDat$diabetes[j]<-diabetes
    pDat$hypertension[j]<-hypertension
    
    rm(hiv,diabetes,hypertension,comorb)
  }
  
  pDat<-pDat %>%
    dplyr::mutate(
      comorb=case_when(
        hiv==1 & diabetes==1 & hypertension==0 ~ "HivDiabetes",
        hiv==1 & diabetes==0 & hypertension==1 ~ "HivHypertension",
        hiv==0 & diabetes==1 & hypertension==1 ~ "DiabetesHypertension",
        hiv==1 & diabetes==1 & hypertension==1 ~ "HivDiabetesHypertension"
      )
    ) %>%
  dplyr::select(!c(rateevent,propHivScreening,propDiabScreening,propHyperScreening))
  
  
}


# function to evaluate ICC
evalIccFromSim <- function(n_cluster, cluster_size, nRegions, nSims, seed) {
  
  set.seed(seed)
  CI_width <- rep(NA, nSims)  # Store CI width from each simulation
  
  for (j in 1:nSims) {
    cat(paste0("Running simulation ", j, " out of ", nSims, "...\n"))
    
    # Simulate data
    new_data <- simulate_data(n_cluster = n_cluster, cluster_size = cluster_size, nRegions = nRegions)
    
    # Fit model
    model <- tryCatch({
      glmer(event_count ~ studyArm + ageCentred + sex + comorb + country + (1 | site), 
            data = new_data, 
            family = poisson(link = "log"))
    }, error = function(e) {
      cat("Error in model fitting:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(model)) next  # Skip this iteration if model fitting fails
    
    # Compute ICC
    icc_results <- tryCatch({
      icc(model, ci = TRUE, iterations = 100)
    }, error = function(e) {
      cat("Error in ICC calculation:", e$message, "\n")
      return(NULL)
    })
    
    # Extract ICC values
    if (!is.null(icc_results) && "ICC_adjusted" %in% names(icc_results) && length(icc_results$ICC_adjusted) >= 3) {
      icc_adj <- icc_results$ICC_adjusted[1]
      icc_LCI <- icc_results$ICC_adjusted[2]  
      icc_HCI <- icc_results$ICC_adjusted[3]  
      CI_width[j] <- icc_HCI - icc_LCI  
    } else {
      CI_width[j] <- NA  
    }
  }
  
  # Compute mean CI width across successful simulations
  mean_CI_width <- mean(CI_width, na.rm = TRUE)
  
  return(mean_CI_width)
}


# calculate precision / CI width
final_results <- evalIccFromSim(n_cluster=4, cluster_size=25, nRegions=3, nSims=100, seed=123)

save(final_results, file = "final_results.rda")
