################################################################################
# Title: Multilevel joint modeling of longitudinal and binary outcomes for 
#        hierarchically structured data
# Copyright: Grace C. Zhou, Seongho Song, Rhonda D. Szczesniak
# Initiation Date: 2021.06.10
# Last Update Date: 2021.11.8
# Notes:
# 1. This is an example code mainly for jitter data set analyzed by spep-JM4 model
# 2. R version for rstan::sampling: 3.6.1(2019-07-05), x86_64-pc-linux-gnu
# 3. R version for the remaining: 4.0.2(2020-06-22), x86_64-w64-mingw32
################################################################################

rm(list=ls()) #Clear Global Environment

#LOAD PACKAGES----

require(dplyr) # For Data management 
require(gridExtra) # For plot side by side
require(ggplot2) # For plot
require(rstan) # For stan  
require(bayesplot) # For sampling diagnostics
require(loo) # For WAIC

#READ DATA----
##Notes
##1. JITTER.DATA: Analysis data
##2. JITTER.TEST: 20% patients from analysis data 
##3. JITTER.TRAIN: First 80% obs. from rest 80% patients of analysis data
##4. JITTER.MASK: Mask remaining 20% obs.  

#Load: JITTER.DATA, JITTER.TEST, JITTER.TRAIN, JITTER.MASK
#Make sure file was saved in your working directory
load('INPUTS.RData') 

##################################################

#DATA DISPLAY (ANALOGOUS FIGURE 1 IN THE ARTICLE)----

AnaDat.sum <- JITTER.DATA %>% 
  group_by(centerID) %>% 
  summarise(ss=n_distinct(eDWID),freq=round(100*sum(pevent)/n(),2))


AnaDat.new <- JITTER.DATA %>% 
  left_join(AnaDat.sum,by='centerID') %>% 
  mutate(centerID.num=centerID,
         centerID.long=paste0('Center ',centerID.num,' : (Sample Size=',ss,')'),
         centerID.bin=paste0('Center ', centerID.num,' : (PE Freq.=',freq,'%)'),
         PE=ifelse(pevent==1,'Yes','No'))
         

set.seed(2021)
subDat<-AnaDat.new %>%
  select(centerID,eDWID) %>%
  group_by(centerID) %>%
  distinct()%>%
  sample_n(3)

subDat1=AnaDat.new %>%
  filter(eDWID %in% subDat$eDWID)

p.long=AnaDat.new %>%
  ggplot(aes(x=timept, y=FEV1)) + 
  labs(x='Time since first PE (yrs)',y='ppFEV1')+
  geom_point(aes(group=eDWID),color='darkgrey') +
  geom_smooth(se=TRUE, method="loess", colour=alpha('blue3',0.7), span=0.5, level = 0.95, size=1) +
  geom_line(data=subDat1, aes(x=timept,y=FEV1,group=eDWID),size=2) +
  facet_wrap(~centerID.long,scales = "fixed",nrow=5,ncol=1)+
  theme_bw() +
  theme(legend.title = element_blank(),legend.position = "none",
        strip.text = element_text(size=20,face="bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))


p.binary=AnaDat.new %>% 
  ggplot(aes(timept,color=PE)) +
  geom_histogram(aes(y=..density..,fill=PE),alpha=0.5,lwd=1.5)+
  geom_density(aes(y=..density..),size=1)+
  scale_fill_brewer(palette='Set2')+
  scale_color_brewer(palette='Set2')+
  labs(x="Time since first PE (yrs)",
       y='Density',color='PE',fill='PE')+
  facet_wrap(~centerID.bin,scales = "fixed",nrow=5,ncol=1)+
  geom_density(alpha=.5, fill="lightgrey") +
  scale_y_continuous(breaks=seq(0,1,by=0.5)) +
  theme_bw()+
  theme(legend.position = c(0.9, 0.08),
        legend.key.size = unit(0.5, "cm"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        strip.text = element_text(size=20,face="bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))


grid.arrange(p.long, p.binary, ncol=2)


#MODEL (ANALOGOUS TABLE 2 & 3 IN THE ARTICLE)----
##Function: Generate a data list prepared for rstan
##Note: Covariates for LME: basefev, BMIpct, F508.HETER
##      Covariates for GLM: timept, F508.HETER, isOnEnzymes.Y

JOINT_DATA<-function(indata,design.X,design.V){
  
  start_pos=pull(indata %>% 
                 group_by(centerID,patID) %>% 
                 tibble::rownames_to_column('start_pos') %>% 
                 slice(1) %>% 
                 ungroup(),start_pos)

  N=length(unique(indata$eDWID))
  Nlev1=length(unique(indata$centerID))
  NperCenter=pull(indata %>%
                  group_by(centerID) %>%
                  summarise(n=n_distinct(eDWID)),n)
  Nobs.perCenter=pull(indata %>%
                      group_by(centerID) %>%
                      tally(),n)
  Nobs.perSub=pull(indata %>%
                   group_by(centerID,patID) %>%
                   tally(),n)
  
  object.x=lm(design.X, indata)
  X=model.matrix(object.x)
  V=model.matrix(design.V, indata)

  joint.data <- list(Nobs = nrow(indata), # total number of obs.
                     N = N, # total number of patients
                     start_pos=c(as.numeric(start_pos),nrow(indata)+1), # starting point of each patient
                     T=Nobs.perSub, # number of obs. per patient
                     t=indata$timept, # time since first PE
                     NpredsX = ncol(X), # number of predictors for LME
                     X = X, # design matrix for LME
                     NpredsV = ncol(V), # number of predictors for GLM
                     V = V, # design matrix for GLM
                     Nlev1=Nlev1, # number of centers (level1)
                     levind1 = rep(1:Nlev1,times=Nobs.perCenter), #id index for centers (level1)
                     levind2 = rep(1:N,times=Nobs.perSub),#id index for patients (level2)
                     y = indata$FEV1, # FEV1 (continuous) outcome
                     r = indata$pevent, # PE (binary) outcome
                     sdscal=sd(residuals(object.x)) #standard deviation
                )
  
  return(joint.data)

}

JD=JOINT_DATA(indata=JITTER.TRAIN, design.X = as.formula(FEV1~1+basefev+BMIpct+F508.HETER),
              design.V = as.formula(FEV1 ~ 1+timept+F508.HETER+isOnEnzymes.Y)) 

##FIT JM4 WITH SPEP
##Warning message:In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :'-E' not found is safe to ignore
Stan.SPEP<-rstan::stan_model(stanc_ret=stanc("SPEP_JM4.stan"),verbose=FALSE) #read in STAN code

##BE AWARE OF 4-5 HRS WAITING TIME
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

FIT.JM4.SPEP<-rstan::sampling(Stan.SPEP,data=JD,
                              warmup=2000,iter=5000,thin=3,chains=2,
                              control=list(adapt_delta=0.95,max_treedepth=15),
                              seed=1781002850, init_r=0.2,
                              save_warmup=FALSE)


##EXTRACT INFO. 
stanfit=FIT.JM4.SPEP
pars=c('alpha','beta','pow_r','sigmalev1','sigmalev2','sigmaeps','rho1','rho2','rho','tau')

print(stanfit,par=pars,prob=c(.025, .975)) # ANALOGOUS TABLE 3 IN THE ARTICLE
rstan::get_elapsed_time(stanfit) # ANALOGOUS TABLE 6 IN WEB APPENDIX I
rstan::get_seed(stanfit) #1781002850 (ANALOGOUS TABLE 6 IN WEB APPENDIX I)

log_lik_y=loo::extract_log_lik(stanfit, parameter_name = "log_lik_y", merge_chains = TRUE) 
log_lik_r=loo::extract_log_lik(stanfit, parameter_name = "log_lik_r", merge_chains = TRUE) 

loo::waic(log_lik_y) #WAIC1=23308.8 (ANALOGOUS TABLE 2 IN THE ARTICLE)
loo::waic(log_lik_r) #WAIC2=1934.0 (ANALOGOUS TABLE 2 IN THE ARTICLE)

##HMC DIAGNOSTICS (ANALOGOUS FIGURE IN WEB APPENDIX G)
rstan::traceplot(stanfit,par=pars) + ggtitle('Traceplot')+theme(plot.title = element_text(hjust = 0.5)) 
rstan::stan_ac(stanfit,pars) + ggtitle('Autocorrelation Plot')+theme(plot.title = element_text(hjust = 0.5)) 
bayesplot::mcmc_neff(neff_ratio(stanfit,pars=pars))+yaxis_text(hjust = 0)+ggtitle('Effective sample size')+theme(plot.title = element_text(hjust = 0.5))

#DIAGNOSTICS (ANALOGOUS FIGURE 2 IN THE SUPPLEMENTARY)----
##Empirical transformed residual(Diggle et al. 2015)

stanfit<-FIT.JM4.SPEP
summary.stanfit=rstan::summary(stanfit)
{
  sigma1.est=summary.stanfit$summary['sigmalev1','mean']
  sigma2.est=summary.stanfit$summary['sigmalev2','mean']
  sigma.est=summary.stanfit$summary['sigmaeps','mean']
  tau.est=summary.stanfit$summary['tau','mean']
  rho.est=summary.stanfit$summary['rho','mean']
}

indata<- OUT.TRAIN
sd.resid.save <- NULL

for(sub in unique(indata$eDWID)){
  
  sub.traind=indata %>%
    filter(eDWID==sub)
  
  time=sub.traind$timept
  raw.resid=sub.traind$FEV1-sub.traind$Xalpha
  
  n_i=nrow(sub.traind)
  J=matrix(1,n_i,n_i)
  I=diag(n_i)
  R=outer(time,time,function(x,y) exp(-abs(x-y)*rho.est))
  
  Vi=(sigma1.est^2+sigma2.est^2) * J + tau.est^2 * R + sigma.est^2 * I
  
  S <- t(chol(Vi))
  sd.resid<- solve(S) %*% raw.resid
  
  sd.resid.save<-c(sd.resid.save,sd.resid)
  
}

##Plot
par(mfrow=c(2,2))

##(a) Standardized residual
plot(indata$Xalpha,sd.resid.save, col="skyblue",xlab="Fitted values",
     ylab="Standardized residuals",type='p',main='')
abline(0,0,col="black")
lines(lowess(indata$Xalpha,sd.resid.save),col='red',lty=3, lwd = 2)

##(b) Standardized residual vs. time 
plot(indata$timept, sd.resid.save, col="skyblue",xlab="Time since baseline(yrs)",
     ylab="Standardized residuals", type="p",)
abline(0,0,col="black")
lines(lowess(indata$timept, sd.resid.save),col="red",lty=3, lwd = 2)

##(c) Kernel density plot
hist(sd.resid.save, freq=F, ylim=c(0,0.5),xlab="Standardized residuals",ylab="Density",main='',col='skyblue')
lines(density(rnorm(30000)),lwd=2)

##(d) Q-Q plot
qqnorm(sd.resid.save,main="",col='skyblue'); qqline(sd.resid.save,lwd=1)

## Add title
mtext("Residual Plots", side = 3, line = -1, outer = TRUE)



