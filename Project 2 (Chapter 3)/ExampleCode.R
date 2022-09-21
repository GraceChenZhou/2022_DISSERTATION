################################################################################
# Title: Multilevel Bayesian joint modeling for a longitudinal marker and time-to-recurrent
# event to characterize heterogeneity in multi-center studies
# Copyright: Grace C. Zhou, Seongho Song, Pedro M. Afonso, Eleni-Rosalina Andrinopoulou, Rhonda D. Szczesniak
# Credit: Samuel Brilleman(Generate standata.list), Pedro M. Afonso(Simulation)
# Date: 2022.06.30
# Notes:
# 1. This is an example code to generate simulated data and fit the data with two-stage method and joint model
# 2. We set association=current value, time scale=calendar
# 3. You can skip two-stage method, but need to pre-calculate y1_eta_q & a_prior_scale
# 4. R version 4.0.5 (2021-03-31)
################################################################################
#Load----
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(survival)
library(lme4)
library(cmdstanr) #v0.4.0
source('ExampleCodeFunctions.R') #some essential functions

#Simulation----

########################################################################
#Function purpose: generate simulated data
#Input# L @ target num center
#Input# n @ target num patient per center
#Input# n_scl @ scale num for patient (in case final sample size is less than target)
#Input# assoc @ choices of association structure: current slope (etaslope) or current value (etavalue)
#Input# time.scale @ choices of risk scale: gap or calendar time

#Output# A data list contains:
#Output# longitudinal data (datLong); event data (datRec);
#Output# seed number (seed);
#Output# association value (assoc);
#Output# true values for longitudinal submodel (Long.info);
#Output# true values for event submodel (Rec.info)
########################################################################
sim_recJData <- function(L = 6, n = 80, n_scl = 1.5, assoc=c('etaslope','etavalue'),
                         time.scale=c('gap','cal')){
  
  seed=sample.int(.Machine$integer.max, 1)
  set.seed(seed)
  N_target <- n * L # target total number of patients
  N <- N_target* n_scl
  tmax <- 10 # maximum follow-up time (type I censoring)
  
  ##############################################################################
  # longitudinal outcome 1/2
  ##############################################################################
  ## parameters true values
  y1_betas <- c("Intercept" = 77, 'time'=-1.27, 'time2'=-0.11)
  y1_aux <- 9.08 # measurement error sd
  u_sd <- c(24.18,3)
  b_sd <- 5.07
  b <- rnorm(L,0,b_sd)
  rho <- 0.158
  S.scl   <- matrix(c(1, rho, rho, 1), nrow=2)
  S   <- diag(u_sd) %*% S.scl %*% diag(u_sd)
  U   <- MASS::mvrnorm(N, mu=c(0,0), Sigma=S)
  ##############################################################################
  # terminal outcome
  ##############################################################################
  ## parameters true values
  if (assoc=='etaslope'){
    
    e_betas <- c("(Intercept)" = -5, "w1" = 0.3, "w2" = 0.2) # phi = exp(Intercept)
    e_aux <- c(1.2, 1.1, 0.92, 0.91, 1.03, 1.04) # shape par. from Weibull
    a_beta <- c(-0.3, -0.28, -0.2, -0.17, -0.13, -0.1)
    
    
  } else if (assoc=='etavalue'){
    
    e_betas <- c("(Intercept)" = -2, "w1" = 0.3, "w2" = 0.2) # phi = exp(Intercept)
    e_aux <- c(1.2, 1.1, 0.92, 0.91, 1.03, 1.04) # shape par. from Weibull
    a_beta <- c(-0.04,-0.05,-0.04, -0.03, -0.02,-0.04)
    
  }
  
  basehaz = rep(e_aux, each=n*n_scl)
  
  w1 <- rbinom(N, 1, 0.5) 
  w2 <- rnorm(N, 0, 1)
  
  W_t <- cbind("(Intercept)" = 1, "w1" = w1, "w2" = w2)
  eta_t <- as.vector(W_t %*% e_betas) 
  
  H <- function(t, u, i) {
    h <- function(s) { 
      log_haz <- log(basehaz[i]) + (basehaz[i] - 1) * log(s) + eta_t[i]
      exp(log_haz) 
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
  }
  
  u_t <- runif(N)
  ter_times <- numeric(N)
  
  for(i in seq_len(N)) {
    root <- try(uniroot(H, interval = c(1E-8, 200),
                        u = u_t[i], i = i)$root, TRUE)  
    ter_times[i] <- if (!inherits(root, "try-error")) root else NA
  }
  
  surv_na <- !is.na(ter_times)
  if(sum(surv_na) < N_target) stop("Terminal: Not enough patients. Increase 'n_scl'.")
  
  datSurv0 <- data.frame(id    = seq_len(N),
                         Ttime = ter_times,
                         w1    = w1,
                         w2    = w2,
                         b0    = rep(b,each=n*n_scl),
                         U0    = U[,1],
                         U1    = U[,2],
                         centerID = rep(seq_len(L),each=n*n_scl),
                         basehaz = basehaz,
                         assoc     = rep(a_beta,each=n*n_scl)) %>% 
    filter(!is.na(Ttime))
  
  
  datSurv <- datSurv0 %>% mutate(Tstatus=as.numeric(Ttime <= tmax),
                                 Stime=pmin(Ttime, tmax))
  
  ##############################################################################
  # recurring outcome
  ##############################################################################
  
  W_r <- model.matrix(~ w1 + w2, datSurv)
  
  v_sd <- 1.5
  frailty=rnorm(nrow(datSurv),0,v_sd)
  
  e_eta <- as.vector(W_r %*% e_betas)
  
  if (assoc=='etaslope'){
    if (time.scale=='gap'){
      H <- function(t, u, i, tstart) {
        h <- function(s) { 
          
          y_eta <- y1_betas[2] + y1_betas[3]*2*(s+tstart) + datSurv$U1[i]
          log_haz <- log(datSurv$basehaz[i]) + (datSurv$basehaz[i]-1) * log(s) + 
            e_eta[i] + y_eta * datSurv$assoc[i] + frailty[i]
          exp(log_haz)
          
        }
        
        integrate(h, lower = 0, upper = t)$value + log(u)
        
      }
    } else if (time.scale=='cal'){
      H <- function(t, u, i, tstart) {
        h <- function(s) { 
          
          y_eta <- y1_betas[2] + y1_betas[3]*2*(s+tstart) + datSurv$U1[i]
          log_haz <- log(datSurv$basehaz[i]) + (datSurv$basehaz[i]-1) * log(s+tstart) + 
            e_eta[i] + y_eta * datSurv$assoc[i] + frailty[i]
          
          exp(log_haz)
          
        }
        
        integrate(h, lower = 0, upper = t)$value + log(u)
        
      }
    } else (stop('Time scales are limited to gap and cal'))
  } else if (assoc=='etavalue'){
    if (time.scale=='gap'){
      H <- function(t, u, i, tstart) {
        h <- function(s) { 
          
          y_eta <- y1_betas[1] + y1_betas[2]*(s+tstart) + y1_betas[3]*(s+tstart)^2 + 
            datSurv$b0[i]+datSurv$U0[i]+datSurv$U1[i]*(s+tstart)
          log_haz <- log(datSurv$basehaz[i]) + (datSurv$basehaz[i]-1) * log(s) + 
            e_eta[i] + y_eta * datSurv$assoc[i] + frailty[i]
          exp(log_haz)
          
        }
        
        integrate(h, lower = 0, upper = t)$value + log(u)
        
      }
    } else if (time.scale=='cal'){
      H <- function(t, u, i, tstart) {
        h <- function(s) { 
          
          y_eta <- y1_betas[1] + y1_betas[2]*(s+tstart) + y1_betas[3]*(s+tstart)^2 + 
            datSurv$b0[i]+datSurv$U0[i]+datSurv$U1[i]*(s+tstart)
          log_haz <- log(datSurv$basehaz[i]) + (datSurv$basehaz[i]-1) * log(s+tstart) + 
            e_eta[i] + y_eta * datSurv$assoc[i] + frailty[i]
          exp(log_haz)
          
        }
        
        integrate(h, lower = 0, upper = t)$value + log(u)
        
      }
    }
  } else (stop('Association structures are limited to etaslope and etavalue'))
  
  stop_times <- start_times <- id_times <- list()
  j <- 1
  for(i in seq_along(datSurv$id)) {
    tstart <- 0
    while(!is.na(tstart) & tstart < datSurv$Stime[i]) {
      u_r <- runif(1)
      root <- try(uniroot(H,interval = c(1E-8, 200),
                          u = u_r, i = i, tstart = tstart)$root, TRUE)  
      delta.t <- if(!inherits(root, "try-error")) root else NA
      start_times[[j]] <- tstart
      stop_times[[j]] <- tstart + delta.t
      dur <- 28/365.25 #later add event duration (random or not-random)
      tstart <- tstart + delta.t + dur
      id_times[[j]] <- datSurv$id[i]
      j <- j + 1
    }
  }
  
  datRec0 <- data.frame(id     = unlist(id_times),                       
                        tstart = unlist(start_times),
                        tstop  = unlist(stop_times))
  
  datRec00 <- datRec0 %>% left_join(datSurv,by='id') %>% 
    mutate(status=as.numeric(!is.na(tstop) & tstop <= Stime),
           tstop=pmin(tstop, Stime, na.rm=TRUE)) 
  
  datRec1 <- datRec00 %>% group_by(id) %>% mutate(ft.tstop=first(tstop)) %>% 
    filter(ft.tstop>1E-3) %>% select(-ft.tstop) %>% ungroup()
  
  print(datRec1 %>% group_by(id) %>% tally() %>% summarize(nmax=max(n)) %>% pull(nmax))
  
  first.choice=unique(datRec1$id)
  
  if (length(first.choice)>=N_target){
    
    datRec <- datRec1 %>% 
      filter(id %in% sample(first.choice,N_target,replace=FALSE)) %>% 
      select(-Ttime,-Tstatus)
    
  } else{
    stop('Rec: Not enough patients. Increase n_scl')
  }
  
  ##############################################################################
  # longitudinal outcome 2/2
  ##############################################################################
  dumLong <- data.frame(id   = unique(datRec$id),
                        time = 0)
  
  datLong0 <- data.frame(id   = datRec$id,
                         time = datRec$tstop)
  
  datLong1 <- rbind(dumLong,datLong0) %>% arrange(id,time)
  
  datRec.info <- datRec %>% select(id,centerID,b0,U0,U1) %>% distinct()
  
  datLong2 <- datLong1 %>% left_join(datRec.info,by='id') %>% arrange(centerID,id,time)
  
  X <- model.matrix(~ time + I(time^2), datLong2)
  
  datLong <- datLong2 %>% 
    mutate(y_eta=X %*% y1_betas + b0 + U0 + U1*time,
           eps=rnorm(nrow(.),0,y1_aux),
           y=y_eta+eps)  %>% 
    select(-b0,-U0,-U1,-eps)
  
  datLong$id  <- match(datLong$id, unique(datLong$id)) # rename IDs
  datRec$id  <- match(datRec$id, unique(datRec$id)) # rename IDs

  ##############################################################################
  # save
  ##############################################################################
  datList=list(datLong = datLong, 
               datRec = datRec, 
               seed=seed,
               assoc=assoc,
               Long.info=c('y1_betas'=y1_betas, 'y1_aux'=y1_aux,'b_sd'=b_sd,
                           'u_sd'=u_sd,'rho'=rho),
               Rec.info=c('e_beta'=e_betas,'e_aux'=e_aux,'a_beta'=a_beta,'v_sd'=v_sd))
  
  remove(datLong0,datLong1,datLong2,datRec.info,datRec0,datRec1,datSurv,datSurv0,dumLong,
         id_times,S, S.scl, start_times,stop_times,U,W_r,W_t,X,b,b_sd,basehaz,delta.t,
         dur,e_aux,e_betas,e_eta,eta_t,first.choice,i,j,root,surv_na,ter_times,
         tstart,u_r,u_sd,u_t,w1,w2,y1_aux,y1_betas)
  
  return(datList)
  
}
########################################################################
simJD=sim_recJData(L = 6, n = 80, n_scl = 1.5, assoc='etavalue', time.scale='cal')
########################################################################
########################################################################
#Function purpose: generate data list for Stan
#Input(default)# formulaLong=y~time+I(time^2)+(1|centerID)+(1+time|id)
#Input(default)# formulaEvent=~w1+w2
#Input# datLong @ longitudinal data
#Input# datRec @ event data
#Input# seed @ seed number
#Input# qnodes @ num nodes (Gauss-Kronrod quadrature)
#Input# assoc @ choices of association structure: current slope (etaslope) or current value (etavalue)
#Input# id_var @ id variable in datRec
#Input# time_var @ time variable in datRec
#Input# time.scale @ choices of risk scale: gap or calendar time
#Output# A data list contains:
#Output# dataQ @ datRec at quadrature points;
#Output# lme.fit @ lme fit results for longitudinal submodel;
#Output# standata.jm @ data list for Stan;
#Output# staninit.jm @ data list for Stan;

get_SIMstandata <- function(datLong, datRec, seed, qnodes=7, assoc=c('etaslope','etavalue'),
                            id_var='id', time_var='time',time.scale=c('gap','cal')){
  
  formulaLong=y~time+I(time^2)+(1|centerID)+(1+time|id)
  formulaEvent=~w1+w2
  min_prior_scale = 1e-12
  
  #standata.long----
  Design.matrices <- make_assoc_parts_for_stan(datLong,formulaLong)
  predictors.X <- Design.matrices$predictors.x
  y1_X.scl <- Design.matrices$x.scl
  y1_Xbar <- Design.matrices$x.bar
  y1_prior_scale0=2.5*sd(datLong$y)
  y1_prior_scale=pmax(min_prior_scale, y1_prior_scale0 / apply(predictors.X, 2L, get_scale_value))
  z <- Design.matrices$z
  
  #standata.surv----
  e_mod <- handle_e_mod(data = datRec, qnodes = qnodes, id_var = id_var)
  datRec.info=datRec %>% select(id,centerID,w1,w2) %>% distinct()
  dataQ <- data.frame(id = e_mod$cids, 
                      time = e_mod$cpts.cal) %>% 
           left_join(datRec.info,by='id')
  
  e_Xq <- model.matrix(formulaEvent, data = dataQ)
  e_Xq <- drop_intercept(e_Xq)
  e_Xbar <- as.vector(colMeans(e_Xq))
  e_Xq.scl <- sweep(e_Xq, 2, e_Xbar, FUN = "-")  # Centre the design matrix
  
  #standata.jm----
  if (assoc=='etavalue'){
    
    y1_Xq_raw <- model.matrix(~time+I(time^2),dataQ)
    y1_Xq_raw <- drop_intercept(y1_Xq_raw)
    
    y1_Xq <- sweep(y1_Xq_raw, 2, y1_Xbar, FUN = "-")
    y1_Zq <- model.matrix(~1+time,dataQ)
    
  } else if (assoc=='etaslope'){
    
    y1_Xq <- model.matrix(~I(time*2),dataQ) 
    y1_Zq <- cbind(0, model.matrix(~1,dataQ))
    
  }
  
  L=length(unique(datRec$centerID))
  
  if (time.scale=='gap') {
    basehaz_X=log(as.vector(e_mod$cpts.gap))
  } else if (time.scale=='cal') {
    basehaz_X=log(as.vector(e_mod$cpts.cal))
  } else 
    stop('Time scales are limited to gap and calendar')
  
  standata=list(
    y_N=nrow(datLong),
    y_K=ncol(predictors.X), # exclude intercept
    y1=as.vector(datLong$y),
    y1_X=y1_X.scl, # exclude intercept+scale
    y1_Xbar=y1_Xbar,
    b_N=L, #number of centers
    b_K=ncol(z$centerID), #intercept
    u_N=as.integer(e_mod$Npat), # number of patients
    u_K=ncol(z$id), #intercept + slope
    y1_C_id=datLong$centerID,
    y1_U_id=datLong$id,
    y1_Z=t(z$id), # b_K times y_N matrix
    e_K=ncol(e_Xq.scl), # x1, x2
    e_N=nrow(datRec),
    a_K=L, #association par. for current value
    Npat=as.integer(e_mod$Npat),
    Nevents=as.integer(e_mod$Nevents),
    qnodes=qnodes,
    Nobs_times_qnodes=nrow(datRec)*qnodes,
    nrow_e_Xq=nrow(dataQ),
    e_Xq=e_Xq.scl, ## exclude intercept Xq.scl
    e_Xbar=e_Xbar,
    basehaz_df=L,
    basehaz_X=basehaz_X,
    norm_const=e_mod$norm_const,
    qwts=as.array(e_mod$qwts),
    y1_Xq=y1_Xq, # exclude intercept
    y1_Zq=t(y1_Zq),
    y1_Cq_id=dataQ$centerID,
    y1_Uq_id=dataQ$id,
    y1_prior_scale=y1_prior_scale, # exclude intercept
    e_prior_scale=rep(10,ncol(e_Xq.scl)), 
    y_prior_scale_for_intercept=100,
    e_prior_scale_for_intercept=20,
    y_prior_scale_for_aux=2.5*sd(datLong$y),
    e_prior_scale_for_aux=rep(5,L),
    b_prior_scale=10,
    u_prior_scale=c(10,10),
    u_prior_df=rep(1,ncol(z$id)),
    u_prior_regularization=2,
    Dev_index=ifelse(assoc=='etaslope',0,1)
  )
  
  #staninit.jm----
  lme.fit<- nlme::lme(lme4::nobars(formulaLong),
                      random = list(centerID=~1,id=~time), #random intercept
                      data = datLong,
                      control = nlme::lmeControl(opt = "optim"),
                      method='ML')
  
  y1_betas <- nlme::fixef(lme.fit)
  VarCor=nlme::VarCorr(lme.fit)
  int4center=which(rownames(VarCor)=='centerID =')+1
  int4id=which(rownames(VarCor)=='id =')+1
  
  b_sd=as.numeric(VarCor[int4center,'StdDev'])
  u_sd=as.numeric(c(VarCor[int4id,'StdDev'],VarCor['time','StdDev']))
  
  y1_aux=as.numeric(VarCor['Residual','StdDev'])
  y1_aux_unscaled=y1_aux/standata$y_prior_scale_for_aux
  
  rho=as.numeric(VarCor['time','Corr'])
  S.scl   <- matrix(c(1, rho, rho, 1), nrow=2)
  S   <- diag(u_sd) %*% S.scl %*% diag(u_sd)
  
  u_cholesky=chol(S)
  
  y1_beta=y1_betas[-1] #drop intercept
  y1_z_beta=y1_beta/standata$y1_prior_scale
  y1_gamma=y1_betas[1]+sum(standata$y1_Xbar * y1_beta)
  
  e_alpha=0.2
  e_beta=rep(0.2,L)
  e_aux=rep(1,L)
  a_beta=rep(0.2,L)
  
  e_gamma=e_alpha-standata$norm_const+sum(standata$e_Xbar,e_beta)
  e_aux_unscaled=e_aux/standata$e_prior_scale_for_aux
  e_z_beta=e_beta/standata$e_prior_scale
  
  set.seed(seed)
  staninit=list(
    y1_gamma=as.vector(y1_gamma),
    e_gamma=e_gamma,
    y1_z_beta=as.vector(y1_z_beta),
    y1_aux_unscaled=y1_aux_unscaled,
    e_aux_unscaled=e_aux_unscaled,
    e_z_beta=e_z_beta,
    a_beta=rep(0.2,L),
    b_sd=b_sd,
    u_sd=u_sd,
    u_cholesky=u_cholesky,
    v_sd=0.2,
    z_b_vec=rnorm(standata$b_N,0,1),
    z_v_vec=rnorm(standata$u_N,0,1),
    z_u_mat=t(MASS::mvrnorm(standata$u_N, mu=c(0,0), Sigma=matrix(c(1,0,0,1),nrow=2,byrow = T))))
  
  #save----
  
  stanList=list(dataQ=dataQ,
                lme.fit=lme.fit,
                standata.jm=standata,
                staninit.jm=staninit)
  
  return(stanList)
  
}

########################################################################
simSD=get_SIMstandata(datLong=simJD$datLong, datRec=simJD$datRec, 
                      seed=simJD$seed, qnodes=7, 
                      assoc='etavalue', time.scale='cal',
                      id_var='id', time_var='time')
########################################################################

#Execute in Stan----

standata.jm <- simSD$standata.jm
staninit.jm<- simSD$staninit.jm
seed <- simJD$seed

##Read in Stan models
###Note: 1. stan file can be edited by .txt or Rstudio
###      2. save files in your working directory
file.long <- file.path("LONG.stan")
long.sim <- cmdstan_model(file.long) # model name

file.event <- file.path("EVENT.stan")
event.sim <- cmdstan_model(file.event) # model name

file.jm <- file.path("JM.stan")
jm.sim <- cmdstan_model(file.jm) # model name

##Two-stage Method

###Stage 1

fit.long <- long.sim$sample(
    data = standata.jm,
    chains = 2, 
    save_warmup = FALSE,
    parallel_chains = 2,
    refresh = 500,
    adapt_delta=0.95,
    max_treedepth=12,
    seed=seed,
    init = function() staninit.jm
  )

###Stage 2
####Note: Need to update y1_eta_q & a_prior_scale from Stage 1 first

  y1_eta_q_extract <- suppressWarnings(fit.long$draws('y1_eta_q',format='df') %>% 
                                       select(-.chain, -.iteration, -.draw))
  standata.jm$y1_eta_q <- unname(colMeans(y1_eta_q_extract))
  
  dataQ <- simSD$dataQ
  dataQ$y1_eta_q=standata.jm$y1_eta_q
  a_prior_scale <- c()
  for (i in unique(dataQ$centerID)){
    subdataQ <- dataQ %>% filter(centerID==i)
    a_beta_scale <- get_scale_value(subdataQ$y1_eta_q)
    a_prior_scale[i] <- pmax(1e-12,2.5/a_beta_scale)
  }
 
  standata.jm$a_prior_scale <- a_prior_scale
  simSD$standata.jm=standata.jm
  staninit.jm$a_z_beta=staninit.jm$a_beta/a_prior_scale
  simSD$staninit.jm=staninit.jm

  fit.event <- event.sim$sample(
      data = standata.jm,
      chains = 2, 
      save_warmup = FALSE,
      parallel_chains = 2,
      refresh = 500,
      adapt_delta=0.95,
      max_treedepth=12,
      seed=seed,
      init = function() staninit.jm
  )
  
##Joint model
###Note: there might be some warnings, but they are safe to ignore as the chain will get converged afterwards  
fit.jm <- jm.sim$sample(
      data = standata.jm, #Note: this is the updated list from Stage 2
      chains = 2, 
      save_warmup = FALSE,
      parallel_chains = 2,
      refresh = 500,
      adapt_delta=0.95,
      max_treedepth=12,
      seed=seed,
      init = function() staninit.jm #Note: this is the updated list from Stage 2
  )

##Read estimates
pars.all <- c('y1_alpha','y1_beta','y1_aux','b_sd','u_sd','rho',
              'e_alpha','e_beta','e_aux','a_beta','v_sd')

fit.jm$summary(c(pars.all)) %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  print(n=nrow(.))

