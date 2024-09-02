# Load required libraries
library(survival)
library(cmprsk)
library(parallel)
library(foreach)

# Set seed for reproducibility
set.seed(35501225)


# Number of patients
n <- 500

# Number of repetitions
num_repetitions <- 100
##do it for 500 replications too

generate_cif_data <- function() {
  
  # Generate a confounder (e.g., L as our measured confounder)
  L <- rnorm(n, mean = 0, sd = 1)
  
  # Generate initial treatment initiation times using a Weibull distribution
  lambda <- 2.5  # Shape parameter for Weibull distribution
  scale_param <- 0.585 * exp(0.2 * L)
  W_star_weibull <- rweibull(n, shape = lambda, scale = scale_param)
  
  summary(W_star_weibull)
  
  # Identify patients whose Weibull-generated treatment times are beyond 0.5
  g <- 0.5
  late_treated_indices <- which(W_star_weibull > g)
  num_7<-round(0.178 * length(late_treated_indices))
  
  sampled_indices <- sample(late_treated_indices, num_7)
  
  W_star_weibull[sampled_indices] <- 7 + runif(num_7, min = 0, max = 3)  # Adding a random time beyond 7 years for variation
  
  
  # Final treatment initiation times
  v <- W_star_weibull
  summary(v)
  hist(v)
  
  # Determine treatment status at any time t
  A <- ifelse(v <= g, 1, 0)
  table(A)
  
  # Time to events
  # Event 1: PCa Death
  #lambda1 <- 0.020
  lambda1 <- 0.030
  betaA1 <- -0.78
  betaL <- 0.20

  
  # Pre-treatment
  eventT10 <- rexp(n, rate = lambda1 * exp(0* L))
  # Post-treatment
  newT1 <- rexp(n, rate = lambda1 * exp(0 * L + betaA1))
  
  # Event time 1
  eT1 <- ifelse(eventT10 < v, eventT10, (newT1 + v))
  
  # Event 2 (not dependent on A_t)
  lambda2 <- 0.037
  betaL2 <- 0.25
  eventT20 <- rexp(n, rate = lambda2 * exp(0 * L))
  eT2 <- eventT20
  
  
  #random censoring
  time_censor <- pmin(7, rexp(n, rate = 0.025)) # Assuming an exponential distribution for censoring times
  summary(time_censor)
  
  
  # Combine both events and censoring times to get the survival time
  followuptime <- pmin(eT1, eT2, time_censor)
  summary(followuptime)
  hist(followuptime)
  
  # Generate event indicator (1: prostate cancer death, 2: other death, 0: censored)
  death <- ifelse(eT1 == followuptime, 1, ifelse(eT2 == followuptime, 2, 0))
  table(death)
  
  # Create a data frame
  Simdata <- data.frame(
    Patientid = 1:n,
    L = L,
    A = A,
    Treatment_time = v,
    Followuptime = followuptime,
    Death = death,
    Vital_status = ifelse(death == 0, 0, 1)
  )
  
  # Convert followuptime and treatment_time to days
  Simdata$Followuptime <- round(Simdata$Followuptime * 365.25, digits = 2)
  Simdata$Treatment_time <- round(Simdata$Treatment_time * 365.25, 2)
  
  # Overall survival at 7 years
  km <- survfit(Surv(Followuptime, Vital_status) ~ A, data = Simdata)
  summary(km, times = 2556)
  
  plot(km, col=c("blue","red"), fun="event", xlab="Time (days)",
       ylim=c(0,0.40), ylab="1-Survival probability",
       lwd=2,conf.int=F,main=c(paste("Complement of the KM Plot - Simulated Data"), 
                               paste("n=500 Patients (Scenario 1) 1 replicate"))
  )
  
  legend(x="topleft",c("Treatment strategy A=0: 34%","Treatment strategy A=1: 28%"),
         col=c("blue","red"),lty=1,lwd=2,bty="n", cex=0.8)
  
  # CIFs at 7 years
  cif <- survfit(Surv(Followuptime, as.factor(Death)) ~ A, data = Simdata)
  summary(cif, times = 2556)
  
  plot(cif,ylim=c(0,0.4),xlab="Time (days)",ylab="Cumulative incidence", 
       col=c("blue","red", "blue", "red"),lwd=2,lty=c(1,1,2,2),conf.int=F,
       main=c(paste("CIF Plot - Simulated Data"), paste("n=500 Patients (Scenario 1) 1 replicate")))
  legend(x="topleft",c("Treated PCa CIF: 8.6%","Control PCa CIF: 11.3%", 
                       "Treated OD CIF: 18.9%", "Control OD CIF: 22.5%"),
         col=c("red","blue", "red", "blue"),lwd=2,lty=c(1,1,2,2),bty="n", c=0.8)
  
  ################################################################################
  ##########CCW TIME
  
  ##Treatment definition: patients should receive radiotherapy within 6 months of diagnosis
  Simdata$rt6<-ifelse(Simdata$Treatment_time<=182.62 & Simdata$A==1, 1, 0)
  table(Simdata$rt6)
  ##260 received treatment within 6 months
  
  
  #STEP 1: Cloning and creation of outcome and follow-up time in each counterfactual arm
  #Lets start by creating clones
  #Clones assigned to the control arm (no radiotherapy within 6 months of diagnosis)
  Simdata_control<-Simdata
  Simdata_control$arm<-"Control"
  
  ##Clones assigned to the treatment arm(received radiotherapy within 6 months of diagnosis)
  Simdata_treatment<-Simdata
  Simdata_treatment$arm<-"treatment"
  
  #STEP 2: Censoring and data manipulation to suit each of the defined treatment strategies
  #for each counterfactual clone
  
  #CONTROL ARM: Patients who received radiotherapy within 6 months in the control arm
  Simdata_control$outcome[Simdata_control$rt6==1]<-0
  Simdata_control$fup[Simdata_control$rt6==1]<-Simdata_control$Treatment_time[Simdata_control$rt6==1]
  
  #CONTROL ARM: patients who do not receive radiotherapy within 6 months but may receive later or not
  Simdata_control$outcome[(Simdata_control$rt6==0 & Simdata_control$Followuptime>182.62) | (Simdata_control$rt6==0 & Simdata_control$Treatment_time>182.62)]<- Simdata_control$Vital_status[(Simdata_control$rt6==0 & Simdata_control$Followuptime>182.62)| (Simdata_control$rt6==0 & Simdata_control$Treatment_time>182.62)]
  Simdata_control$fup[(Simdata_control$rt6==0 & Simdata_control$Followuptime>182.62) | (Simdata_control$rt6==0 & Simdata_control$Treatment_time>182.62)]<- Simdata_control$Followuptime[(Simdata_control$rt6==0 & Simdata_control$Followuptime>182.62) | (Simdata_control$rt6==0 & Simdata_control$Treatment_time>182.62)]
  
  ###CONTROL ARM: Patients who die or lost to follow up within 6 months of diagnosis without having radiotherapy
  Simdata_control$outcome[Simdata_control$rt6==0 & Simdata_control$Followuptime<=182.62]<-Simdata_control$Vital_status[Simdata_control$rt6==0 & Simdata_control$Followuptime<=182.62]
  Simdata_control$fup[Simdata_control$rt6==0 & Simdata_control$Followuptime<=182.62]<-Simdata_control$Followuptime[Simdata_control$rt6==0 & Simdata_control$Followuptime<=182.62]
  
  
  ##TREATMENT ARM: Patients received radiotherapy within 6 months of diagnosis
  Simdata_treatment$outcome[Simdata_treatment$rt6==1]<-Simdata_treatment$Vital_status[Simdata_treatment$rt6==1]
  Simdata_treatment$fup[Simdata_treatment$rt6==1]<-Simdata_treatment$Followuptime[Simdata_treatment$rt6==1]
  
  ##TREATMENT ARM: patients who do not receive radiotherapy but are still alive or at risk
  Simdata_treatment$outcome[(Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime>182.62) | (Simdata_treatment$rt6==0 & Simdata_treatment$Treatment_time>182.62)]<-0
  Simdata_treatment$fup[(Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime>182.62) | (Simdata_treatment$rt6==0 & Simdata_treatment$Treatment_time>182.62)]<-182.62
  
  ###TREATMENT ARM: Patients who die or lost to follow up within 6 months of diagnosis without having radiotherapy
  Simdata_treatment$outcome[Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime<=182.62]<-Simdata_treatment$Vital_status[Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime<=182.62]
  Simdata_treatment$fup[Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime<=182.62]<-Simdata_treatment$Followuptime[Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime<=182.62]
  
  #Censoring indicators (artificial censoring)
  ###CONTROL ARM CENSORING: Patients who received radiotherapy within 6 months in the control arm
  Simdata_control$censoring[Simdata_control$rt6==1]<-1
  
  ####CONTROL ARM CENSORING: patients who do not receive radiotherapy within 6 months but may receive later or not
  Simdata_control$censoring[Simdata_control$rt6==0 | (Simdata_control$rt6==0 & Simdata_control$Treatment_time>182.62)]<- 0
  
  ##CONTROL ARM CENSORING: Patients who died or were lost to follow up within 6 months and did not receive radiotherapy
  Simdata_control$censoring[Simdata_control$rt6==0 & Simdata_control$Followuptime<=182.62]<-0
  
  
  ###TREATMENT ARM CENSORING
  ##Patients who receive treatment within 6 months of diagnosis
  Simdata_treatment$censoring[Simdata_treatment$rt6==1]<-0
  
  ##TREATMENT ARM: patients who do not receive radiotherapy but are still alive or at risk
  Simdata_treatment$censoring[(Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime>182.62) | (Simdata_treatment$rt6==0 & Simdata_treatment$Treatment_time>182.62)]<-1
  
  ###TREATMENT ARM: Patients who die or lost to follow up within 6 m of diagnosis without having radiotherapy
  Simdata_treatment$censoring[Simdata_treatment$rt6==0 & Simdata_treatment$Followuptime<=182.62]<-0
  
  ##Bind both the datasets for the control and treated clones
  Allsimdata<-rbind(Simdata_control, Simdata_treatment)
  
  #STEP 3: Weighting step
  # we will begin by reshaping the data to a long format then perform the weighting at each step
  # The control clones will be weighted using a Cox model
  # The treated clones will be weighted using a logistic regression model
  
  ##Long format
  t_events<-sort(unique(Allsimdata$fup))
  times<-data.frame("tevent"=t_events, "timeID_t"=seq(1:length(t_events)))
  
  ####Let's start with the treatment clone (radiotherapy within 6 MONTHS of diagnosis)
  ##Create an entry variable  for everyone
  Simdata_treatment$Tstart<-0
  
  ##Split the dataset at each followuptime of event until the event happens and sort it
  data.long_t6<-survSplit(Simdata_treatment, cut=t_events, end="fup", start="Tstart", event="outcome", id="timeID")
  data.long_t6<-data.long_t6[order(data.long_t6$timeID, data.long_t6$fup),]
  
  #####Splitting the original dataset at each followuptime of censoring and sorting it
  data.long_t.cens6<- survSplit(Simdata_treatment, cut=t_events, end="fup", start="Tstart", event = "censoring", id="timeID")
  data.long_t.cens6<-data.long_t.cens6[order(data.long_t.cens6$timeID, data.long_t.cens6$fup),]
  
  ###Replace the censoring variable in data.long_t with that from data.long_t.cens
  data.long_t6$censoring<-data.long_t.cens6$censoring
  
  remove(data.long_t.cens6)
  
  ###Create T-stop
  data.long_t6$Tstop<-data.long_t6$fup
  
  ##Merge and sort
  data.long_t6<-merge(data.long_t6, times, by.x="Tstart", by.y = "tevent", all.x = T)
  data.long_t6<-data.long_t6[order(data.long_t6$timeID, data.long_t6$fup),]
  data.long_t6$timeID_t[is.na(data.long_t6$timeID_t)]<-0
  
  ####LET'S CREATE A COLUMN FOR PROSTATE DEATH AND OTHER DEATH
  data.long_t6$prostateDeath<- data.long_t6$outcome==1 & data.long_t6$Death==1
  data.long_t6$otherDeath<- data.long_t6$outcome==1 & data.long_t6$Death==2
  
  table(data.long_t6$prostateDeath)
  table(data.long_t6$otherDeath)
  
  ##Control clones
  ############################LET'S REPEAT THE SAME FOR THE CONTROL ARM (No radiotherapy group)
  ##Create an entry variable  for everyone
  Simdata_control$Tstart<-0
  
  ##Split the dataset at each followuptime of event until the event happens and sort it
  data.long_c6<-survSplit(Simdata_control, cut=t_events, end="fup", start="Tstart", event="outcome", id="timeID")
  data.long_c6<-data.long_c6[order(data.long_c6$timeID, data.long_c6$fup),]
  
  #####Splitting the original dataset at each followuptime of censoring and sorting it
  data.long_c.cens6<- survSplit(Simdata_control, cut=t_events, end="fup", start="Tstart", event = "censoring", id="timeID")
  data.long_c.cens6<-data.long_c.cens6[order(data.long_c.cens6$timeID, data.long_c.cens6$fup),]
  
  ###Replace the censoring variable in data.long_t with that from data.long_t.cens
  data.long_c6$censoring<-data.long_c.cens6$censoring
  
  remove(data.long_c.cens6)
  
  ###Create T-stop
  data.long_c6$Tstop<-data.long_c6$fup
  
  ##Merge and sort
  data.long_c6<-merge(data.long_c6, times, by.x="Tstart", by.y = "tevent", all.x = T)
  data.long_c6<-data.long_c6[order(data.long_c6$timeID, data.long_c6$fup),]
  data.long_c6$timeID_t[is.na(data.long_c6$timeID_t)]<-0
  
  ##LET'S CREATE A COLUMN FOR PROSTATE DEATH AND OTHER DEATH
  data.long_c6$prostateDeath<- data.long_c6$outcome==1 & data.long_c6$Death==1
  data.long_c6$otherDeath<- data.long_c6$outcome==1 & data.long_c6$Death==2
  
  table(data.long_c6$prostateDeath)
  table(data.long_c6$otherDeath)
  
  ###FINAL DATASET
  data_final<-rbind(data.long_t6, data.long_c6)
  data_final<-merge(data_final, times, by="timeID_t", all.x = T)
  data_final<-data_final[order(data_final$timeID, data_final$fup),]
  
  ##tevent has NA's which we can replace with 0's for data_final
  data_final$tevent[is.na(data_final$tevent)]<-0
  data_final<-data_final[order(data_final$Followuptime),]
  
  data_final$compet<-ifelse(data_final$outcome==1 & data_final$Death==1, 1, 0)
  data_final$compet<-ifelse(data_final$outcome==1 & data_final$Death==2, 
                            2, data_final$compet)
  table(data_final$compet)
  
  
  ########################################
  #####################################
  #################################### WITH L ONLY
  ####Weighting for the treated clones
  data.long_t6<-data_final[data_final$arm=="treatment",]
  ##Let's fit the censoring model
  newdatat<-with(Simdata_treatment, data.frame(censoring, L,
                                               fup, outcome, Treatment_time, 
                                               Followuptime, arm, rt6,Patientid))
  ##LOGISTIC REGRESSION MODEL FITTED FOR EACH ROW OF DATA: 1 PATIENT PER ROW
  
  ms_cens_t<-glm(censoring~ L , data=Simdata_treatment, family = binomial)
  summary(ms_cens_t) 
  
  newdatat$prob<-predict(ms_cens_t,newdata=newdatat, type="response")
  newdatat$probunc<-1-newdatat$prob
  
  
  #Numerator for stabilised weights
  ms_cens_tnum<-glm(censoring~ 1 , data=Simdata_treatment, family = binomial)
  summary(ms_cens_tnum) 
  
  newdatat$probnum<-predict(ms_cens_tnum,newdata=newdatat, type="response")
  newdatat$probuncnum<-1-newdatat$probnum
  
  
  
  
  newdatat$weight<- newdatat$probuncnum / newdatat$probunc
  summary(newdatat$weight)
  
  ##Time to merge so that weights before 6 months are set to 1
  str(newdatat)
  newdatat<-newdatat[,-1:-8]
  newdatat<-newdatat[,-2:-5]
  data.long_t6<-merge(data.long_t6, newdatat)
  data.long_t6$weight[data.long_t6$fup<=182.62]<-1
  
  summary(data.long_t6$weight)
  
  data.long_t6<-data.long_t6[order(data.long_t6$Patientid, data.long_t6$fup),]
  
  data.long_t6<-data.long_t6[order(data.long_t6$timeID, data.long_t6$fup),]
  
  
  wtrt<-survfit(Surv(Tstart, Tstop, as.factor(compet)) ~ arm, 
                data=data.long_t6, id=data.long_t6$Patientid, weights = data.long_t6$weight)
  summary(wtrt, times=2556)
  
  
  
  ###Control clones dataset
  ###LET'S GO TO THE CONTROL ARM
  data.long_c6<-data_final[data_final$arm=="Control",]
  
 
  ##Cox model
  ms_cens_c6<-coxph(Surv(Tstart, Tstop, censoring)~L, ties = "efron",
                    data= data.long_c6)
  summary(ms_cens_c6)
  
  # Check proportional hazards assumption using Schoenfeld residuals
  cox.zph(ms_cens_c6)
  plot(cox.zph(ms_cens_c6))
  
  ######How I've done it
  hazard<-basehaz(ms_cens_c6, centered = F)
  names(hazard)<-c("hazard", "t")
  haz0<-data.frame(0,0)
  names(haz0)<-c("hazard", "t")
  hazard<-rbind(hazard, haz0)
  hazard<-unique(merge(hazard, times, by.x = "t", by.y = "tevent", all.x = T))
  hazard$timeID_t[is.na(hazard$timeID_t)]<-0
  data.long_c6<-merge(data.long_c6, hazard, by="timeID_t",all.x = T)
  beta<-coef(ms_cens_c6)
  data.long_c6$lin_pred<- as.matrix(data.long_c6[,c(4)])%*%beta
  
  data.long_c6$probunc<-exp(-(data.long_c6$hazard)*exp(data.long_c6$lin_pred))
  summary(data.long_c6$probunc)
  
  #Numerator for stabilised weights
  ms_cens_c6num<-coxph(Surv(Tstart, Tstop, censoring)~1, ties = "efron",
                       data= data.long_c6)
  summary(ms_cens_c6num)
  ######How I've done it
  haznum<-basehaz(ms_cens_c6num, centered = F)
  names(haznum)<-c("hazard", "t")
  haz0<-data.frame(0,0)
  names(haz0)<-c("hazard", "t")
  haznum<-rbind(haznum, haz0)
  haznum<-unique(merge(haznum, times, by.x = "t", by.y = "tevent", all.x = T))
  haznum$timeID_t[is.na(haznum$timeID_t)]<-0
  data.long_c6<-merge(data.long_c6, haznum, by="timeID_t",all.x = T)
  betanum<-coef(ms_cens_c6num)
  
  data.long_c6$probuncnum<-exp(-(data.long_c6$hazard.y)*exp(0))
  summary(data.long_c6$probuncnum)
  
  
  
  
  data.long_c6$weight<- data.long_c6$probuncnum/data.long_c6$probunc
  summary(data.long_c6$weight)
  
  
  #data.long_c6<-data.long_c6[order(data.long_c6$Patientid, data.long_c6$fup),]
  #data.long_c6<-data.long_c6[order(data.long_c6$timeID, data.long_c6$fup),]
  
  
  
  
  wctrl<-survfit(Surv(Tstart, Tstop, as.factor(compet)) ~ arm, 
                 data=data.long_c6, id=data.long_c6$Patientid, weights = data.long_c6$weight)
  summary(wctrl, times=2556)
  
  
  
  
  ###LET'S PUT ALL THE CIF'S INTO ONE PLOT
  summary_cif<-summary(wctrl)
  summary_cif
  # Store the CIFs in a data frame
  cif_data_cont <- data.frame(Time = summary_cif$time, CIF = summary_cif$pstate, 
                              SE=summary_cif$std.err)
  
  ##Treated
  summary_cif2<-summary(wtrt)
  # Extract the cumulative incidence probabilities
  summary_cif2
  # Store the CIFs in a data frame
  cif_data_treat <- data.frame(Time = summary_cif2$time, CIF = summary_cif2$pstate)
  
  
  library(ggplot2)
  # Assuming you have calculated CIFs for treatment groups and stored them as data frames, 
  #like cif_treatment1 and cif_treatment2
  
  cif_data_cont$Time_years <- cif_data_cont$Time / 365.25
  cif_data_treat$Time_years <- cif_data_treat$Time / 365.25
  
  cif_data_cont$CIF.1_percent <- cif_data_cont$CIF.1 * 100
  cif_data_cont$CIF.2_percent <- cif_data_cont$CIF.2 * 100
  cif_data_treat$CIF.1_percent <- cif_data_treat$CIF.1 * 100
  cif_data_treat$CIF.2_percent <- cif_data_treat$CIF.2 * 100
  
  
  cif.func.cont1<-stepfun(cif_data_cont$Time_years,c(0,cif_data_cont$CIF.1_percent))
  
  cif.func.cont1(seq(0,7,0.1))
  
  # Create step functions for CIFs
  cont_stepfun1_adj <- stepfun(cif_data_cont$Time_years,c(0,cif_data_cont$CIF.1_percent))
  cont_stepfun1_adj(seq(0,7,0.1))
  cont_stepfun2_adj <- stepfun(cif_data_cont$Time_years,c(0,cif_data_cont$CIF.2_percent))
  cont_stepfun2_adj(seq(0,7,0.1))
  
  treat_stepfun1_adj <- stepfun(cif_data_treat$Time_years, c(0, cif_data_treat$CIF.1_percent))
  treat_stepfun1_adj(seq(0,7,0.1))
  treat_stepfun2_adj <- stepfun(cif_data_treat$Time_years, c(0, cif_data_treat$CIF.2_percent))
  treat_stepfun2_adj(seq(0,7,0.1))
  
  
  
  
  
  
  ##################################################
  #################################################
  ########################################WITHOUT L
  
  ###LET'S PUT ALL THE CIF'S INTO ONE PLOT
  summary_cif_a<-summary(cif)
  summary_cif_a
  # Store the CIFs in a data frame
  cif_data_original <- data.frame(Time = summary_cif_a$time, CIF = summary_cif_a$pstate, Treatment= summary_cif_a$strata)
  
  cif_data_cont_a<-subset(cif_data_original, cif_data_original$Treatment=="A=0")
  cif_data_treat_a <-subset(cif_data_original, cif_data_original$Treatment=="A=1")
  
  
  
  cif_data_cont_a$Time_years_a <- cif_data_cont_a$Time / 365.25
  cif_data_treat_a$Time_years_a <- cif_data_treat_a$Time / 365.25
  
  cif_data_cont_a$CIF.1_percent_a <- cif_data_cont_a$CIF.1 * 100
  cif_data_cont_a$CIF.2_percent_a <- cif_data_cont_a$CIF.2 * 100
  cif_data_treat_a$CIF.1_percent_a <- cif_data_treat_a$CIF.1 * 100
  cif_data_treat_a$CIF.2_percent_a <- cif_data_treat_a$CIF.2 * 100
  
  treat_stepfun1_unadj <- stepfun(cif_data_treat_a$Time_years_a, c(0,cif_data_treat_a$CIF.1_percent_a))
  treat_stepfun1_unadj(seq(0,7,0.1))
  cont_stepfun1_unadj <- stepfun(cif_data_cont_a$Time_years_a, c(0,cif_data_cont_a$CIF.1_percent_a))
  cont_stepfun1_unadj(seq(0,7,0.1))
  treat_stepfun2_unadj <- stepfun(cif_data_treat_a$Time_years_a, c(0,cif_data_treat_a$CIF.2_percent_a))
  treat_stepfun2_unadj(seq(0,7,0.1))
  cont_stepfun2_unadj <- stepfun(cif_data_cont_a$Time_years_a, c(0,cif_data_cont_a$CIF.2_percent_a))
  cont_stepfun2_unadj(seq(0,7,0.1))
  
  
  
  
  return(list(treat_stepfun1_adj(seq(0,7,0.1)), cont_stepfun1_adj(seq(0,7,0.1)), 
              treat_stepfun2_adj(seq(0,7,0.1)), cont_stepfun2_adj(seq(0,7,0.1)),
              treat_stepfun1_unadj(seq(0,7,0.1)), cont_stepfun1_unadj(seq(0,7,0.1)), 
              treat_stepfun2_unadj(seq(0,7,0.1)),
              cont_stepfun2_unadj(seq(0,7,0.1))))
  #,              Simdata_replications=sim_reps))
  
  
  #}
  
  
}




###CHECK FOR BIAS, MSE, COVERAGE
#Create an empty list to store Simdata
#sim_reps <- list() 

# Repeat the data generation process for the desired number of times
cif_data_list <- replicate(num_repetitions, generate_cif_data(), simplify = FALSE)


# Initialize an empty list to store data frames for the cifs
cif_dfs <- list()


# Access CIF step functions for each repetition
for (i in 1:num_repetitions) {
  cat("Repetition:", i, "\n")
  cif_stepfun_list <- cif_data_list[[i]]
  # Check cif_stepfun_list[[1]], ...
}


# Extract CIFs from all repetitions for each index from 1 to 8
for (i in 1:8) {
  cif_data <- sapply(cif_data_list, "[[", i)
  cif_dfs[[i]] <- data.frame(cif_data)
  
  ##colnames
  col_names <- paste("CIF", i, rep(seq(1:num_repetitions), each = 1), sep = "_")
  colnames(cif_dfs[[i]]) <- col_names
}


# Combine the data frames into one data frame
final_df <- do.call(cbind, cif_dfs)

##Simdata replicates
#Simdata_replications <- as.data.frame(cif_stepfun_list$Simdata_replications)
#Simlist<-cif_stepfun_list[[9]]





### For the average rowmeans we have
# Initialize an empty list to store data frames
average_cif_dfs <- list()

# Calculate row-wise means for each CIF
for (i in 1:8) {
  # Extract CIFs for the current CIF number
  cif_cols <- grep(paste("CIF", i, sep = "_"), colnames(final_df))
  cif_data <- final_df[, cif_cols]
  
  # Calculate row-wise means
  average_cif_dfs[[i]] <- rowMeans(cif_data, na.rm = TRUE)
}

# Combine the data frames into one data frame
average_df <- do.call(cbind, average_cif_dfs)

# Set column names for the average data frame
colnames(average_df) <- paste("Average_CIF", 1:8, sep = "_")
average_df <- as.data.frame(average_df)

# Save the data frame as a CSV file
write.csv(final_df, "final_cif_data_scenario1.csv", row.names = FALSE)
write.csv(average_df, "average_cif_data_scenario1.csv", row.names = FALSE)

##LET'S LOAD THE DATA FOR PLOTTING

#import the average_cif_data_scenario1.csv
#average_df<-average_cif_data_scenario1
#final_df<-final_cif_data_scenario1

##Plot the cif's for nsim=100 and n=500
library(dbplyr)
library(tidyverse)

average_df <- rename(average_df, "Average PCa treated*" = Average_CIF_1,
                     "Average PCa control*"=Average_CIF_2,
                     "Average OD treated*"=Average_CIF_3,
                     "Average OD control*"=Average_CIF_4,
                     "Average PCa treated"=Average_CIF_5,
                     "Average PCa control"=Average_CIF_6,
                     "Average OD treated"=Average_CIF_7,
                     "Average OD control"=Average_CIF_8)
average_df$Time<-seq(0,7,0.1)
df_long <- gather(average_df, key = "Variables", value = "Values", -Time)

# Create separate dataframes for the two subsets of variables
df_subset1 <- df_long[df_long$Variables %in% c("Average PCa treated*", "Average PCa control*",
                                               "Average OD treated*", "Average OD control*"), ]
df_subset2 <- df_long[df_long$Variables %in% c("Average PCa treated", "Average PCa control",
                                               "Average OD treated", "Average OD control"), ]


if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)
# Plot variables against time for the second subset
last_p1 <- df_subset1 %>%
  group_by(Variables) %>%
  filter(Time == max(Time)) %>%
  ungroup()  

# Determine the nudge direction based on the difference
#between the maximum and minimum values
last_p1<- last_p1 %>%
  mutate(nudge_y = ifelse(max(Values) > min(Values), 1, -1))


p1 <- ggplot(df_subset1, aes(x = Time, y = Values, linetype = Variables, color = Variables)) +
  geom_line() +
  labs(x = "Time in Years", y = "CIF Estimates in %") +
  scale_color_manual(values = c("blue", "red", "blue", "red")) +  # Use blue and red colors
  scale_linetype_manual(values = c("dashed", "dashed","solid", "solid")) +  # Use dashed and solid lines
  theme_minimal() + # Use a minimal theme
  scale_y_continuous(limits = c(0, 40)) +  # Set limits for the y-axis
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  labs(color = "Variables", title = "Risk Estimates Controlling for L Using CCW Approach \nScenario 1: n=500 patients and 100 replicates")+
  geom_text_repel(data = last_p1, aes(label = round(Values, digits=1)), vjust = ifelse(last_p1$nudge_y == -1, -0.01, 0.501), direction = "y", size = 2, show.legend = FALSE)  # Add text labels with exact values at the end, suppressing the legend

#geom_text(data = last_p1, aes(label = round(Values, 1)), vjust = 0, hjust = -0.1, size = 2, show.legend = FALSE)  # Add text labels with exact values at the end


p1




# Plot variables against time for the second subset
last_p2 <- df_subset2 %>%
  group_by(Variables) %>%
  filter(Time == max(Time)) %>%
  ungroup()  

# Determine the nudge direction based on the difference
#between the maximum and minimum values
last_p2 <- last_p2 %>%
  mutate(nudge_y = ifelse(max(Values) > min(Values), 1, -1))



p2 <- ggplot(df_subset2, aes(x = Time, y = Values, linetype = Variables, color = Variables)) +
  geom_line() +
  labs(x = "Time in Years", y = "CIF Estimates in %") +
  scale_color_manual(values = c("blue", "red", "blue", "red")) +  # Use blue and red colors
  scale_linetype_manual(values = c("dashed", "dashed","solid", "solid")) +  # Use dashed and solid lines
  theme_minimal() + # Use a minimal theme
  scale_y_continuous(limits = c(0, 40)) +  # Set limits for the y-axis
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  labs(color = "Variables", 
       title = "Risk Estimates Ignoring L Using Standard Approach \nScenario 1: n=500 patients and 100 replicates")+
  geom_text_repel(data = last_p2,aes(label = sprintf("%.1f", Values)), vjust = ifelse(last_p2$nudge_y == -1, -0.01, 0.501), direction = "y", size = 2, show.legend = FALSE)  # Add text labels with exact values at the end, suppressing the legend

p2





p1+theme(legend.position = "bottom", legend.text = element_text(size = 4.6)) 

p2+theme(legend.position = "bottom", legend.text = element_text(size = 4.6)) 


##BIAS
##bias


s1<-summary(cif1)
truecif1<- data.frame(Time = s1$time, CIF = s1$pstate)
truecif1$Time<-round(truecif1$Time/365.25,digits=2)
truecif1<-truecif1[,-c(2)]
truecif1$CIF.1<-truecif1$CIF.1*100
truecif1$CIF.2<-truecif1$CIF.2*100
truecif1 <- rename(truecif1, "True PCa treated" = CIF.1,
                   "True OD treated"=CIF.2)
# Create step functions for CIFs
treat_truePCA<- stepfun(truecif1$Time, c(0, truecif1$"True PCa treated"))
True_PCa_treated<-treat_truePCA(seq(0,7,0.1))
treat_trueOD <- stepfun(truecif1$Time, c(0, truecif1$"True OD treated"))
True_OD_treated<-treat_trueOD(seq(0,7,0.1))

truetreat<-data.frame(True_PCa_treated, True_OD_treated)

newc1<- gather(truecif1, key = "Variables", value = "Values", -Time)

s2<-summary(cif2)
truecif2<- data.frame(Time = s2$time, CIF = s2$pstate)
truecif2$Time<-round(truecif2$Time/365.25,digits=2)
truecif2<-truecif2[,-c(2)]
truecif2$CIF.1<-truecif2$CIF.1*100
truecif2$CIF.2<-truecif2$CIF.2*100
truecif2 <- rename(truecif2, "True PCa control" = CIF.1,
                   "True OD control"=CIF.2)
# Create step functions for CIFs
control_truePCA<- stepfun(truecif2$Time, c(0, truecif2$"True PCa control"))
True_PCa_control<-control_truePCA(seq(0,7,0.1))
control_trueOD <- stepfun(truecif2$Time, c(0, truecif2$"True OD control"))
True_OD_control<-control_trueOD(seq(0,7,0.1))

truecontrol<-data.frame(True_PCa_control, True_OD_control)
newc2<- gather(truecif2, key = "Variables", value = "Values", -Time)

allc<-rbind(newc1, newc2)

allaverage<-cbind(average_df, truetreat, truecontrol)

allaverage$bias1_t_CCW<-allaverage$"Average PCa treated*" - allaverage$True_PCa_treated
allaverage$bias1_c_CCW<-allaverage$"Average PCa control*" - allaverage$True_PCa_control
allaverage$bias2_t_CCW<-allaverage$"Average OD treated*" - allaverage$True_OD_treated
allaverage$bias2_c_CCW<-allaverage$"Average OD control*" - allaverage$True_OD_control
allaverage$bias1_t_STD<-allaverage$"Average PCa treated" - allaverage$True_PCa_treated
allaverage$bias1_c_STD<-allaverage$"Average PCa control" - allaverage$True_PCa_control
allaverage$bias2_t_STD<-allaverage$"Average OD treated" - allaverage$True_OD_treated
allaverage$bias2_c_STD<-allaverage$"Average OD control" - allaverage$True_OD_control

a2<-allaverage[c(1:9)]



# Reshape the data
long_data <- a2 %>%
  pivot_longer(cols = -Time, 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(method = if_else(str_detect(variable, "\\*$"), "CCW", "Standard")) %>%
  mutate(variable = str_remove(variable, "\\*$")) %>%
  rename(Time = Time, A = variable, Estimates_CIF = value, Method = method)

a3<-allaverage[c(9:13)]

a3<- a3 %>%
  rename("Average PCa treated"=True_PCa_treated, "Average OD treated"=True_OD_treated,
         "Average PCa control"=True_PCa_control, "Average OD control"=True_OD_control)

long_data2 <- a3 %>%
  pivot_longer(cols = -Time, 
               names_to = "variable", 
               values_to = "value") %>%
  rename(Time = Time, A = variable, True_CIF = value)


dataforbias<-merge(long_data,long_data2, by=c("Time", "A"))

dataforbias$bias<-dataforbias$Estimates_CIF - dataforbias$True_CIF

# Plot bias over time
ggplot(dataforbias, aes(x = Time, y = bias, color = Method, group = Method))+
  geom_line() +
  geom_point(size=0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add a horizontal line at y = 0
  facet_grid(Method~A)+
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  scale_y_continuous(breaks = seq(-9, 9, by = 3), limits = c(-9, 9)) +  # Enforce y-axis limits
  ggtitle("Bias Over Time for the CCW and the Standard Approach \nScenario 1: n=500 Patients with 100 Replicates") +
  labs(x = "Time (Years)", y = "Bias in %") +
  theme_minimal()+
  theme(strip.text = element_text(size = 8.5))+
  theme(legend.position = "bottom", legend.text = element_text(size = 5.6))




###THE GRAPH BELOW IS TO PLOT THE TRUTH AND THE CCW AND STD APPROACHES
cols <- c( "blue", "red", "blue","red", "skyblue", "black", "skyblue", "black")

comdata_ccw<-rbind(df_subset1, allc)
ggplot(comdata_ccw, aes(x=Time, y=Values, col=Variables, linetype = Variables)) + 
  geom_line()+
  scale_colour_manual(values = cols)+
  scale_linetype_manual(values = c("dashed", "dashed","solid", "solid","dashed", "dashed","solid", "solid")) +  # Use dashed and solid lines
  theme_minimal() + # Use a minimal theme
  scale_y_continuous(limits = c(0, 40)) +  # Set limits for the y-axis
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  labs(color = "Variables", 
       title = "Comparison of the True Risk Estimates to those Controlling \nfor L Using CCW Approach Scenario 1: n=500 patients \nwith 100 replicates")+
  theme(legend.position = "right", legend.text = element_text(size =4.6)) 


comdata_std<-rbind(df_subset2, allc)
ggplot(comdata_std, aes(x=Time, y=Values, col=Variables, linetype=Variables)) + 
  geom_line()+
  scale_colour_manual(values = cols)+
  scale_linetype_manual(values = c("dashed", "dashed","solid", "solid","dashed", "dashed","solid", "solid")) +  # Use dashed and solid lines
  theme_minimal() + # Use a minimal theme
  scale_y_continuous(limits = c(0, 40)) +  # Set limits for the y-axis
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  labs(color = "Variables", 
       title = "Comparison of the True Risk Estimates to those Ignoring L \nUsing Standard Approach Scenario 1: n=500 patients \nwith 100 replicates")+
  theme(legend.position = "right", legend.text = element_text(size = 4.6)) 



summary_df <- dataforbias %>%
  group_by(Method,A) %>%
  summarise(
    count = n(),                  
    mean_B = mean(bias, na.rm = TRUE),   
    med_B<-median(bias, na.rm=TRUE),
    sd_B<-sd(bias, na.rm=TRUE), 
    min_B = min(bias, na.rm = TRUE),
    max_B = max(bias, na.rm = TRUE)     
    
  )
print(summary_df)






