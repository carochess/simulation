# Load required libraries
library(survival)

# Seed
set.seed(605501225)

# Number of patients
n <- 100000

# Confounder (L as our measured confounder)
L <- rnorm(n, mean = 0, sd = 1)

# Treatment initiation times using a Weibull distribution
lambda <- 2.5  # Shape parameter for Weibull distribution
scale_param <- 0.585 * exp(0.2 * L)
W_star_weibull <- rweibull(n, shape = lambda, scale = scale_param)

hist(W_star_weibull)
summary(W_star_weibull)


##THE ABOVE CODE RETURNS TIMES LESS THAN 0.5 FOR ROUGHLY 51% OF THE PATIENTS

g <- 0.5


# Sample treatment times until all are greater than 0.5
for (i in 1:n) {
  while (TRUE) {
    time <- rweibull(1, shape = lambda, scale = scale_param[i])
    if (time > 0.5) {
      W_star_weibull[i] <- time
      break
    }
  }
}



#the ones who were never treated at all
num_7<-round(0.178 * length(W_star_weibull))
sampled_indices <- sample(length(W_star_weibull), num_7)
W_star_weibull[sampled_indices]  <- 7 + runif(num_7, min = 0, max = 3)  # Adding a random time beyond 7 years for variation



W_star<-W_star_weibull

hist(W_star)


##the above is dependent on L

##
v<-sample(W_star, n, replace = T)
summary(v)
hist(v)

# Determine treatment status at any time t
A <- ifelse(v <= g, 1, 0)
table(A)



##Time to events
#Event 1: PCa Death
lambda1<-0.030; betaA1 <- -0.78; betaL <- 0.95
#pre-treatment
eventT10<-rexp(n, rate = lambda1 * exp(0*L))
hist(eventT10)
#post-treatment
newT1<-rexp(n, rate = lambda1 * exp(0*L + betaA1))
#event time 1
eT1<-ifelse(eventT10<v, eventT10, (newT1+v))
#hist(eT1)
#eT1<-ifelse(eventT10<v, eventT10, (eventT10+v))
#eT1<-eventT10
hist(eT1)


##Event 2 (not dependent on A_t)
lambda2<-0.037; betaL2 <- 1.15
eventT20<-rexp(n, rate = lambda2 * exp(0*L))
hist(eventT20)

eT2<-eventT20

#eT2<-ifelse(eventT20<v, eventT20, (eventT20+v ))
#hist(eT2)

# survival time
followuptime <- pmin(eT1, eT2,7)
hist(followuptime)

# Generate event indicator (1: prostate cancer death, 2: other death, 0: censored)
death <- ifelse(eT1 == followuptime, 1, ifelse(eT2 == followuptime, 2,0))
table(death)

# Create a data frame
Simdata2 <- data.frame(Patientid = 1:n,
                       L=L,
                       A=A,
                       Treatment_time=v,
                       Followuptime=followuptime,
                       Death=death,
                       Vital_status = ifelse(death==0, 0, 1)
)

head(Simdata2)

##Let's convert followuptime and treatment_time to days 
Simdata2$Followuptime<-round(Simdata2$Followuptime*365.25, digits = 2)
Simdata2$Treatment_time<-round(Simdata2$Treatment_time*365.25, 2)

##Overall survival at 7 years
km2 <- survfit(Surv(Followuptime, Vital_status) ~ 0, data = Simdata2)
summary(km2, times=2556)

##CIFs at 7 years
cif2 <- survfit(Surv(Followuptime, as.factor(Death)) ~ 0, data = Simdata2)
summary(cif2, times=2556)


plot(cif2,ylim=c(0,0.40),xlab="Time (days)",ylab="Cumulative incidence", 
     col=c("blue","blue"),lwd=2,lty=c(1,2),conf.int=F,
     main=c(paste("True CIF Plot - Simulated Data Controls"), paste("n=100,000 Patients (Scenario 1)")))
legend(x="topleft",c("Control PCa CIF: 15.9%", "Control OD CIF: 25.6%"),
       col=c("blue", "blue"),lwd=2,lty=c(1,2),bty="n", c=0.8)


###combined true graphs
plot(cif1,ylim=c(0,0.40),xlab="Time (days)",ylab="Cumulative incidence", 
     col=c("red","red"),lwd=2,lty=c(1,2),conf.int=F,
     main=c(paste("True CIF Plot - Simulated Data"), paste("n=100,000 Patients (Scenario 1)")))
lines(cif2,col=c("blue","blue"),lwd=2,lty=c(1,2),conf.int=F)

legend(x="topleft",c("Treated PCa CIF: 8.6%","Control PCa CIF: 10.3%",
                     "Treated OD CIF: 21.7%", "Control OD CIF: 21.3%"),
       col=c("red","blue", "red", "blue"),lwd=2,lty=c(1,1,2,2),bty="n", c=0.8)

