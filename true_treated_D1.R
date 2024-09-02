# Load required libraries
library(survival)

# Seed
set.seed(585501225)

# Number of patients
n <- 100000

# Confounder (L as our measured confounder)
L <- rnorm(n, mean = 0, sd = 1)

# Treatment initiation times using a Weibull distribution
lambda <- 2.5  # Shape parameter for Weibull distribution
scale_param <- 0.585 * exp(0.2 * L)
W_star_weibull <- rweibull(n, shape = lambda, scale = scale_param)

summary(W_star_weibull)

##THE ABOVE CODE RETURNS TIMES LESS THAN 0.5 FOR ROUGHLY 51% OF THE PATIENTS


# Identify patients whose Weibull-generated treatment times are beyond 0.5
g <- 0.5
late_treated_indices <- which(W_star_weibull > g)
num_late <- length(late_treated_indices)

num_7<-round(0.178 * length(late_treated_indices))
sampled_indices <- sample(late_treated_indices, num_7)
W_star_weibull[sampled_indices] <- 7 + runif(num_7, min = 0, max = 3)  # Adding a random time beyond 7 years for variation

# Function to generate valid treatment initiation times 
for (i in 1:n) {
  while (TRUE) {
    time <- rweibull(1, shape = lambda, scale = scale_param[i])
    if (time <= g) {
      W_star_weibull[i] <- time
      break
    }
  }
}

# Generate treatment initiation times

W_star <- W_star_weibull

# Sample treatment initiation times randomly
v <- sample(W_star, n, replace = TRUE)

# Determine treatment status at any time t
A<-ifelse(v<=g, 1, 0)
table(A) ##all patients were untreated
hist(v) ##treatment was within 6 months


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
hist(eT1)

##Event 2 (not dependent on A_t)
lambda2<-0.037; betaL2 <- 1.15
eventT20<-rexp(n, rate = lambda2 * exp(0*L))
hist(eventT20)

eT2<-eventT20


# survival time
followuptime <- pmin(eT1, eT2,7)
hist(followuptime)

# Generate event indicator (1: prostate cancer death, 2: other death, 0: censored)
death <- ifelse(eT1 == followuptime, 1, ifelse(eT2 == followuptime, 2,0))
table(death)

# Create a data frame
Simdata1 <- data.frame(Patientid = 1:n,
                       L=L,
                       A=A,
                       Treatment_time=v,
                       Followuptime=followuptime,
                       Death=death,
                       Vital_status = ifelse(death==0, 0, 1)
)

head(Simdata1)

##Let's convert followuptime and treatment_time to days 
Simdata1$Followuptime<-round(Simdata1$Followuptime*365.25, digits = 2)
Simdata1$Treatment_time<-round(Simdata1$Treatment_time*365.25, 2)

##Overall survival at 7 years
km1 <- survfit(Surv(Followuptime, Vital_status) ~ 1, data = Simdata1)
summary(km1, times=2556)

##CIFs at 7 years
cif1 <- survfit(Surv(Followuptime, as.factor(Death)) ~ 1, data = Simdata1)
summary(cif1, times=2556)

plot(cif1,ylim=c(0,0.40),xlab="Time (days)",ylab="Cumulative incidence", 
     col=c("red","red"),lwd=2,lty=c(1,2),conf.int=F,
     main=c(paste("True CIF Plot - Simulated Data Treated"), paste("n=100,000 Patients (Scenario 1)")))
legend(x="topleft",c("Treated PCa CIF: 8.6%", "Treated OD CIF: 21.7%"),
       col=c("red", "red"),lwd=2,lty=c(1,2),bty="n", c=0.8)


