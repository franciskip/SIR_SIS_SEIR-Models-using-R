library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
setwd("/Users/fyego/Downloads/fk_staff/JUPYTER R")
getwd()

dataprov3 <- read.table("https://bit.ly/2vDqAYN", header = TRUE)
summary(dataprov3)

dataprov3$time <- dataprov3$day
dataprov3$number_infected <- dataprov3$cases



dataprov3 %>% 
  select(time, number_infected) %>% 
  ggplot(aes(x=time, y=number_infected)) +

  geom_line(linetype="dotdash", color="darkorchid") +  geom_point(color="darkorchid", shape=5, show.legend = FALSE) +
  xlab("Time") +
  ylab("Number infected, provided") +
  labs(
    subtitle= "Epidemic Curve, Modelling") +
  theme(legend.position="bottom")

# SIR MODEL FUNCTION 
sir_model <- function(time, state, parameters) { 
    with(as.list(c(state, parameters)), {
      N <- S + I + R
      lambda <- beta * I/N
      dS <- -(lambda * S)
      dI <- (lambda * S) -(gamma * I)
      dR <- (gamma * I)
      return(list(c(dS, dI, dR)))
    })
}

# The sum-of-squares value (SSQ)  FUNCTION (assumes sir_model model already loaded)
# parameters - vector with named elements for sir_model
# idat - df or list containing vectors number_infected (I) and time
SIR_SSQv2 <- function(func, parameters, idat) {  
# MODEL OUTPUT (solving the differential equations):
oderesult <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = func,
                            parms = parameters))
# sum-of-squares (SSQ) of model fit requires idat
idat <- na.omit(idat)
deltas2 <- (oderesult$I[oderesult$time %in% idat$time] - idat$number_infected)^2
SSQ   <- sum(deltas2)
return(SSQ)
}

# initial_state_values and times
initial_state_values <- c(S = 499, I = 1, R = 0)
times <- seq(from = 0, to = 20, by = 1) 

# choose values to start your optimisation
beta_start  <- 0.9
gamma_start <- 0.1

# run optim
(optimised <- optim(par = c(beta = beta_start
                  , gamma = gamma_start)
                  , fn  = SIR_SSQv2
                  , func = sir_model
                  , idat = dataprov3
  ))

# examine optim() output and plot "best" model against example dataset
opt_mod <- as.data.frame(ode(y = initial_state_values  # named vector of initial state values
                  , times = times         # vector of times
                  , func = sir_model      # your predefined SIR function
                  , parms = optimised$par
  ))

nicesubtitle <- "Automated Least-Squares Calibration: Using optim() c2w3_3"

# Compare model I to infection numbers provided
g=opt_mod %>% 
  left_join(dataprov3) %>%
  select(time, I, number_infected) %>% 
  melt(id = "time") %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line(linetype="dotdash") +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values = c("red","darkorchid")) + 
  scale_shape_manual(values = c(4,5)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title=paste(" beta of", round(optimised$par[1],3),
    " gamma of", round(optimised$par[2],3),
    " SSQ", round(optimised$value,0),
    " R0 of", round(optimised$par[1]/optimised$par[2],3)
    ),
    color = "Compartment",
    subtitle= nicesubtitle) +
  theme(legend.position="bottom")  
g

nicesubtitle <- "SIR Model v1a9 SIR Basic Model, check pvacc = pc of .92 \n all-or-nothing vaccine with 75% efficacy"


# MODEL INPUTS:
N        <- 500       # population size
duration <- 14       # total number of days
tsteps   <- 1         # chunk in days  
beta     <- 1.88188642103645       # infection rate day^-1
gamma    <- 0.29465259706606      # recovery rate day^-1
R0 <- beta / gamma


(parameters <- c(
  beta = beta,          # infection rate
  gamma = gamma,        # recovery rate for untreated I
  R0 = R0
  ))  

1-(1/R0)

0.843427002941134/.75

# initial_state_values and times
initial_state_values <- c(S = (N-1) * .88, # 8% are susceptible
                          I = 1, 
                          R = (N-1) * .12) # 92% are vaccinated / immune
# TIMESTEPS:
times <- seq(from = 0, to = duration, by = tsteps)

# SIR MODEL FUNCTION 
sir_model <- function(time, state, parameters) { 
    with(as.list(c(state, parameters)), {
      N <- S + I + R
      lambda <- beta * I/N
      dS <- -(lambda * S)
      dI <- (lambda * S) -(gamma * I)
      dR <- (gamma * I)
      return(list(c(dS, dI, dR)))
    })
}

# MODEL OUTPUT (solving the differential equations):
output2 <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters)) %>% 
  mutate(still_Su = round(S/(S + I + R),digits=5)*100,
    preval_Inf = round(I/(S+ I + R),digits=5)*100,
    propor_Re = round(R/(S + I + R),digits=5)*100,
    Reff = (parameters['beta']/parameters['gamma']) 
    * (S/(S + I + R)) ) 
  

print("peak infection day when I is at its max: ")
output2 %>%
  arrange(-I, time) %>%
  head(1)

print("point when R eff goes below 1: ")
output2 %>%
  filter(Reff < 1) %>% 
  arrange(time) %>%
  head(1)

print("last record for the run: ")
output2 %>%
  arrange(time) %>%
  tail(1)

print("Plotting the number of people in each compartment over time")
output2 %>% 
  filter(time <= 5) %>% 
  select(-still_Su, -preval_Inf, -propor_Re, -Reff) %>% 
  melt(id = "time")  %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line() +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values = c("blue","red","green")) + 
  scale_shape_manual(values = c(0,4,1)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title=paste(" beta of", parameters['beta'],
    " gamma of", parameters['gamma'],
    " R0 of", round(parameters['R0'],3)),
    color = "Compartment",
    subtitle= nicesubtitle) +
  theme(legend.position="bottom") 

print("Plotting the proportion of people in each compartment over time")
## [1] "Plotting the proportion of people in each compartment over time"
output2 %>% 
  filter(time <= 5) %>% 
  select(-still_Su, -preval_Inf, -propor_Re, -Reff) %>% 
  melt(id = "time")  %>% 
  mutate(proportion = value / sum(initial_state_values)) %>% 
  ggplot(aes(x=time, y=proportion, color=variable, shape=variable)) +
  geom_line() +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values = c("blue","red","green")) + 
  scale_shape_manual(values = c(0,4,1)) +
  xlab("Time (days)") +
  ylab("Proportion of people") +
  labs(title=paste(" beta of", parameters['beta'],
    " gamma of", parameters['gamma'],
    " R0 of", round(parameters['R0'],3)),
    color = "Compartment",
    subtitle= nicesubtitle) +
  theme(legend.position="bottom")

output2 %>% 
  filter(time <= 5) %>% 
  select(time, Reff) %>% 
  ggplot(aes(x=time, y=Reff)) +
  geom_line(linetype="dotdash", color="orange") +
  geom_jitter(color="orange", shape=6, show.legend = FALSE) +
  xlab("Time (days)") +
  ylab("R Effective") +
  labs(title=paste(" beta of", parameters['beta'],
    " gamma of", parameters['gamma'],
    " R0 of", round(parameters['R0'],3)),
    color = "Compartment",
    subtitle= nicesubtitle) +
  theme(legend.position="bottom")

# 1 of 4
# Load the flu dataset of reported cases
dataprov4 <- read.csv("SIR/flu_reported.csv")

dataprov4 %>% 
  select(time, number_reported) %>% 
  ggplot(aes(x=time, y=number_reported)) +
  geom_line(linetype="dotdash", color="cadetblue4") +
  geom_point(color="chartreuse3", shape=10, show.legend = FALSE) +
  xlab("Time (days)") +
  ylab("Number of Reported Cases") +
  labs(
    subtitle= "Calibrating SIR model using likelihood c2w4_1") +
  theme(legend.position="bottom")

nicesubtitle <- "SIR Model v1a10 SIR Basic Model, reported data with provided beta and gamma c2w4_1"

print("initial state values and parameters")
# MODEL INPUTS:
N        <- 763      # population size
duration <- 14       # total number of days
tsteps   <- 0.1      # chunk in days  
beta     <- 1.7      # infection rate day^-1
gamma    <- 0.45     # recovery rate day^-1
R0 <- beta / gamma

(parameters <- c(
  beta = beta,          # infection rate
  gamma = gamma,        # recovery rate
  R0 = R0
  ))  

# initial_state_values and times
initial_state_values <- c(S = (N-1),
                          I = 1, 
                          R = 0)

# TIMESTEPS:
times <- seq(from = 0, to = duration, by = tsteps)

# SIR MODEL FUNCTION 
sir_model <- function(time, state, parameters) { 
    with(as.list(c(state, parameters)), {
      N <- S + I + R
      lambda <- beta * I/N
      dS <- -(lambda * S)
      dI <- (lambda * S) -(gamma * I)
      dR <- (gamma * I)
      return(list(c(dS, dI, dR)))
    })
}

# MODEL OUTPUT (solving the differential equations):
output2 <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters)) %>% 
  mutate(still_Su = round(S/(S + I + R),digits=5)*100,
    preval_Inf = round(I/(S+ I + R),digits=5)*100,
    propor_Re = round(R/(S + I + R),digits=5)*100,
    Reff = (parameters['beta']/parameters['gamma']) 
    * (S/(S + I + R)) ) 

output2 %>% 
  left_join(dataprov4) %>%
  select(time, I, number_reported) %>% 
  melt(id = "time") %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line(linetype="dotdash") +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values = c("red","cadetblue4")) + 
  scale_shape_manual(values = c(4,10)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title=paste(" beta of", parameters['beta'],
    " gamma of", parameters['gamma'],
    " R0 of", round(parameters['R0'],3)
    ),
    color = "Compartment",
    subtitle= nicesubtitle) +
  theme(legend.position="bottom")



# Code block 4:
(LL <- sum(dpois(x =  dataprov4$number_reported, lambda = 0.6 * output2$I[output2$time %in%  dataprov4$time],
                 log = TRUE)))

plot(dpois(x =  dataprov4$number_reported, lambda = 0.6 * output2$I[output2$time %in%  dataprov4$time], log = TRUE))

# 2 of 4
# DISTANCE FUNCTION
loglik_function <- function(parameters, idat) { # param  values and dataset
   beta <- parameters[1]  # first value in params is beta
   gamma <- parameters[2] # second value in params is gamma
    
# Simulate the model with initial conditions and timesteps
oderesult2 <- as.data.frame(ode(y = initial_state_values,
  times = times,
  func = sir_model,
  parms = c(beta = beta, # beta fr params of loglik_func
    gamma = gamma)))  # gamma from parameters of loglik_func
    
# Calculate log-likelihood, accounting for the reporting rate of 60%:
LL <- sum(dpois(x = idat$number_reported, lambda = 0.6 * oderesult2$I[oderesult2$time %in% idat$time], log = TRUE))
   return(LL) 
}

# 3 of 4
# OPTIMISATION:
# optim(par = c(parameters['beta'],parameters['gamma']), # doesn't quite get exact same answer
# optim(par = c(1.7, 0.45),
# optim(par = c(1.0, 0.01),
(optimised2 <- optim(par = c(1.7, 0.1), # starting values for beta and gamma - you should get the same result no matter what but not true
      fn = loglik_function,
      idat = dataprov4,
      control = list(fnscale=-1))  # look for max
)

nicesubtitle <- "Automated Maximum Likelihood Calibration: Using optim() c2w4_2"

# Load the flu dataset of the total number infected
dataprov5 <- read.csv("SIR/flu_data1.csv")

# plot total number infected versus prior reported infected
dataprov5 %>% 
  left_join(dataprov4) %>%
  select(time, number_infected, number_reported) %>% 
  melt(id = "time") %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line(linetype="dotdash",show.legend = FALSE) +
  geom_point() +
  scale_color_manual(values = c("cadetblue4","darkorchid")) + 
  scale_shape_manual(values = c(10,5)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(subtitle= "Reported vs. Full Dataset of Infected - c2w4_2") +
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())

# simulate model with the best-fitting (max-likelihood) parameters from prior steps
(parameters <- c(beta = optimised2$par[1],
  gamma = optimised2$par[2],
  R0 = optimised2$par[1]/optimised2$par[2]
  ))

# MODEL OUTPUT (solving the differential equations):
output2 <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters)) %>% 
  mutate(still_Su = round(S/(S + I + R),digits=5)*100,
    preval_Inf = round(I/(S+ I + R),digits=5)*100,
    propor_Re = round(R/(S + I + R),digits=5)*100,
    Reff = (parameters['beta']/parameters['gamma']) 
    * (S/(S + I + R)) ) 

output2 %>% 
  left_join(dataprov5) %>%
  left_join(dataprov4) %>%
  select(time, I, number_infected, number_reported) %>% 
  melt(id = "time") %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line(linetype="dotdash",show.legend = FALSE) +
  geom_point() +
  scale_color_manual(values = c("red", "cadetblue4", "darkorchid")) + 
  scale_shape_manual(values = c(4, 10, 5)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title= paste0("Model *I*nfected vs. Reported vs. Full Dataset of Infected",
    "\n beta of ", round(parameters['beta'],3),
    " gamma of ", round(parameters['gamma'],3),
    " R0 of ", round(parameters['R0'],3)),
    subtitle= nicesubtitle) +
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())  

# MODEL OUTPUT (solving the differential equations):
output2 <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters)) %>% 
  mutate(still_Su = round(S/(S + I + R),digits=5)*100,
    preval_Inf = round(I/(S+ I + R),digits=5)*100,
    propor_Re = round(R/(S + I + R),digits=5)*100,
    Reff = (parameters['beta']/parameters['gamma']) 
    * (S/(S + I + R)) ) 

l=output2 %>% 
  left_join(dataprov4) %>%
  select(time, I, number_reported) %>% 
  melt(id = "time") %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line(linetype="dotdash",show.legend = FALSE) +
  geom_point() +
  scale_color_manual(values = c("red", "cadetblue4", "darkorchid")) + 
  scale_shape_manual(values = c(4, 10, 5)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title= paste0("Model *I*nfected vs. Reported vs. Full Dataset of Infected",
    "\n beta of ", round(parameters['beta'],3),
    " gamma of ", round(parameters['gamma'],3),
    " R0 of ", round(parameters['R0'],3)),
    subtitle= nicesubtitle) +
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())
l

# MODEL OUTPUT (solving the differential equations):
output2 <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters)) %>% 
  mutate(still_Su = round(S/(S + I + R),digits=5)*100,
    preval_Inf = round(I/(S+ I + R),digits=5)*100,
    propor_Re = round(R/(S + I + R),digits=5)*100,
    Reff = (parameters['beta']/parameters['gamma']) 
    * (S/(S + I + R)) ) 

output2 %>% 
  left_join(dataprov5) %>%
  left_join(dataprov4) %>%
  select(time, I, number_infected, number_reported) %>% 
  melt(id = "time") %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line(linetype="dotdash",show.legend = FALSE) +
  geom_point() +
  scale_color_manual(values = c("red", "cadetblue4", "darkorchid")) + 
  scale_shape_manual(values = c(4, 10, 5)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title= paste0("Model *I*nfected vs. Reported vs. Full Dataset of Infected",
    "\n beta of ", round(parameters['beta'],3),
    " gamma of ", round(parameters['gamma'],3),
    " R0 of ", round(parameters['R0'],3)),
    subtitle= nicesubtitle) +
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())  

z=output2 %>% 
  left_join(dataprov4) %>%
  select(time, I, number_reported) %>% 
  melt(id = "time") %>%
  filter(!is.na(value)) %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line(linetype="dotdash") +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values = c("red","cadetblue4")) + 
  scale_shape_manual(values = c(4,10)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title=paste(" beta of", parameters['beta'],
    " gamma of", parameters['gamma'],
    " R0 of", round(parameters['R0'],3)
    ),
    color = "Compartment",
    subtitle= nicesubtitle) +
  theme(legend.position="bottom")
z

library(gridExtra)
options(warn=-1)
grid.arrange(g, z, nrow = 1)


