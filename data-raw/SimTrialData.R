###############
#load packages#
###############

library(tidyverse)
library(lme4)
library(VGAM)
library(Hmisc)
library(gam)

s  <- gam::s
s_ <- VGAM::s 
select <- dplyr::select 
fill <- VGAM::fill

##############
#read in data#
##############

path.data <- "C:\\Users\\John\\Documents\\Work\\Computing\\SASFiles\\SASUniversityEdition\\myfolders\\data\\CATIEDATA2\\CATIEDATA\\"
catie <- read_csv(paste(path.data,"catie_wgt.csv",sep="")) %>% 
  group_by(CATIEID) %>%
  arrange(CATIEID,time) %>%
  mutate(treat.grp=factor(treat.grp,levels=1:5),
         age.grp=factor(age.grp,levels=1:5),
         race=factor(race,levels=1:3),
         phase.change.vis=ifelse(time %in% c(0,1),0,
                             ifelse(PHASE!=lag(PHASE),1,0))) %>%
  ungroup()

#baseline records
catie.0 <- catie %>% filter(time==0) 

##############################
# STEP 0: id and time vectors#
##############################

#get vectors
CATIEID.vec <- unique(catie$CATIEID)
time.vec    <- unique(catie$time)

#simulate person-time observations
n <- length(CATIEID.vec)
obs <- length(time.vec)

CATIEID <- rep(CATIEID.vec,obs)
time <- sort(rep(time.vec,n))

##################################
# STEP 1: simulate baseline covs #
##################################

td           <- rep(rbinom(n,1,mean(catie.0$td)),obs)
zprcort      <- rep(rbinom(n,1,mean(catie.0$zprcort)),obs)
pmat.race    <- matrix(rep(table(catie.0$race)/nrow(catie.0),n*obs),ncol=length(table(catie.0$race)),nrow=n,byrow=TRUE)
race         <- factor(rep(rMultinom(pmat.race,1),obs),levels=1:3)
pmat.age.grp <- matrix(rep(table(catie.0$race)/nrow(catie.0),n*obs),ncol=length(table(catie.0$race)),nrow=n,byrow=TRUE)
age.grp      <- factor(rep(rMultinom(pmat.race,1),obs),levels=1:5)
educ.bin     <- rep(rbinom(n,1,mean(catie.0$educ.bin)),obs)
unemployed   <- rep(rbinom(n,1,mean(catie.0$unemployed)),obs)
exacer       <- rep(rbinom(n,1,mean(catie.0$exacer)),obs)
site.ro      <- rep(rbinom(n,1,mean(catie.0$site.ro)),obs)
site.sh      <- rep(rbinom(n,1,mean(catie.0$site.sh)),obs)
site.uc      <- rep(rbinom(n,1,mean(catie.0$site.uc)),obs)
site.va      <- rep(rbinom(n,1,mean(catie.0$site.va)),obs)

step1.trial <- data.frame(CATIEID,time,td,zprcort,race,age.grp,educ.bin,unemployed,exacer,site.ro,site.sh,site.uc,site.va)
step1.trial.0 <- step1.trial[which(step1.trial$time==0),]

##############################
# STEP 2: simulate treatment #
##############################

fit.treat.grp  <- vgam(formula=treat.grp~zprcort+td,data=catie.0,family=multinomial())
pmat.treat.grp <- predict.vgam(fit.treat.grp,newdata=step1.trial.0,type="response") 
treat.grp      <- factor(rep(rMultinom(pmat.treat.grp,1),obs),levels=1:5) 

step2.trial <- data.frame(step1.trial,treat.grp)

#########################################
# STEP 3: simulate symptom trajectories #
#########################################

fit.cs14 <- lmer(cs14~time+(time|CATIEID)+(time|treat.grp),catie)
cs14 <- predict(fit.cs14,newdata=step2.trial,type="response")

fit.cs16 <- lmer(cs16~time+(time|CATIEID)+(time|treat.grp),catie)
cs16 <- predict(fit.cs16,newdata=step2.trial,type="response")

fit.calg1 <- lmer(calg1~time+(time|CATIEID)+(time|treat.grp),catie)
calg1 <- predict(fit.calg1,newdata=step2.trial,type="response")

fit.weight <- lmer(weight~time+(time|CATIEID)+(time|treat.grp),catie)
weight <- predict(fit.weight,newdata=step2.trial,type="response")

fit.epsmean <- lmer(epsmean~time+(time|CATIEID)+(time|treat.grp),catie)
epsmean <- predict(fit.epsmean,newdata=step2.trial,type="response")

fit.qoltot <- lmer(qoltot~time+(time|CATIEID)+(time|treat.grp),catie)
qoltot <- predict(fit.qoltot,newdata=step2.trial,type="response")

step3.trial <- data.frame(step2.trial,cs14,cs16,calg1,weight,epsmean,qoltot)

#########################################
# STEP 4: simulate outcome trajecotries #
#########################################

fit.pansstotal <- lmer(pansstotal~time+(time|CATIEID)+(time|treat.grp)+cs14+cs16+weight+epsmean+qoltot,catie)
pansstotal <- predict(fit.pansstotal,newdata=step3.trial,type="response")

step4.trial <- data.frame(step3.trial,pansstotal)

#############################################
# STEP 5: simulate phase (treatment) change #
#############################################

#simulate separately for each time and treatment arm

pred.phase.change.vis <- function (fit.data,pred.data,...) {
  
  group_var <- quos(...)
  
  vars <- c(names(pred.data),"phase.change.vis","source")
  pred.data$source <- "pred"
  pred.data$phase.change.vis <-NA
  fit.data$source <- "fit"
  
  fit.data <- fit.data[,vars]
  pred.data <- pred.data[,vars]
  
  fxn.data <- bind_rows(fit.data,pred.data)
  
  mod.pred <- function (indata) {
    
    indata.fit <- indata %>% filter(source=="fit")
    
    #fit <- gam(data=indata.fit,
    #  phase.change.vis ~ 
    #    site.ro + site.sh + site.uc + site.va
    #  + exacer + race + unemployed + age.grp
    #  + calg1 + cs14 + epsmean + pansstotal + weight,
    #  family=binomial("logit")
    #  )

    fit <- gam(data=indata.fit,
      phase.change.vis ~ 
        s(time,4)
      + site.ro + site.sh + site.uc + site.va
      + exacer + race + unemployed + age.grp
      + s(calg1,4) + s(cs14,4) + s(epsmean,4) + s(pansstotal,4) + s(weight,4)
      + s(calg1*time,4) + s(cs14*time,4) + s(epsmean*time,4) + s(pansstotal*time,4) + s(weight*time,4),
      family=binomial("logit")
      )

    print(fit)

    indata.pred <- indata %>% filter(source=="pred")

    X1 <- predict(fit,newdata=indata.pred,type="response")
    X0 <- 1-X1
    pmat <- data.frame(X1,X0)
    names(pmat) <- c("1","0")
    output <- rMultinom(pmat,1) %>% data.frame
    
  }
  
  result <- fxn.data %>% group_by(!!! group_var) %>% do(mod.pred(.)) %>% ungroup() 

}

temp <- pred.phase.change.vis(catie,step4.trial,treat.grp)
colnames(temp) <- c("treat.grp","phase.change.vis")

step5.trial <- data.frame(step4.trial,temp$phase.change.vis) %>% rename(phase.change.vis=temp.phase.change.vis)

#########################################
# STEP 6: create variables for analysis #
#########################################

step6.trial <- step5.trial %>%
  group_by(CATIEID) %>%
  arrange(CATIEID,time) %>%
  mutate(white=ifelse(race==1,1,0),
         black=ifelse(race==2,1,0),
         other=ifelse(race==3,1,0),
         age.grp.1824=ifelse(.data$age.grp==1,1,0),
         age.grp.2534=ifelse(.data$age.grp==2,1,0),
         age.grp.3544=ifelse(.data$age.grp==3,1,0),
         age.grp.4554=ifelse(.data$age.grp==4,1,0),
         age.grp.5567=ifelse(.data$age.grp==5,1,0),
         Bpansstotal=first(.data$pansstotal),
         Bcs14=first(.data$cs14),
         Bcs16=first(.data$cs16),
         Bcalg1=first(.data$calg1),
         Bqoltot=first(.data$qoltot),
         Bepsmean=first(.data$epsmean),
         Chg.pansstotal=.data$pansstotal-first(.data$pansstotal),
         pct.gain=10*(weight-first(.data$weight))/first(.data$weight),
         phase.change.vis=ifelse(time %in% c(0,1),0,
                                 ifelse(phase.change.vis=="1",1,
                                        ifelse(phase.change.vis=="0",0,NA))),
         phase.change.cum=cumsum(.data$phase.change.vis),
         phase.change.cum.rec=ifelse(.data$phase.change.cum>=1,
                                     1,
                                     .data$phase.change.cum),
         lead.pansstotal=lead(.data$pansstotal),
         treat.grp.ola=ifelse(.data$treat.grp==1,1,0),
         treat.grp.que=ifelse(.data$treat.grp==2,1,0),
         treat.grp.ris=ifelse(.data$treat.grp==3,1,0),
         treat.grp.per=ifelse(.data$treat.grp==4,1,0),
         treat.grp.zip=ifelse(.data$treat.grp==5,1,0)
         ) %>%
  ungroup()

############################
# STEP 7: simulate dropout #
############################

#simulate separately for each time and treatment arm

pred.studydisc <- function (fit.data,pred.data,...) {
  
  group_var <- quos(...)
  
  vars <- c(names(pred.data),"studydisc","source")
  pred.data$source <- "pred"
  pred.data$studydisc <-NA
  fit.data$source <- "fit"
  
  fit.data <- fit.data[,vars]
  pred.data <- pred.data[,vars]
  
  fxn.data <- bind_rows(fit.data,pred.data)
  
  mod.pred <- function (indata) {
    
    indata.fit <- indata %>% filter(source=="fit")
    
    #fit <- gam(data=indata.fit,
    #           studydisc ~ 
    #             Bpansstotal + site.ro + site.sh + site.uc + site.va
    #           + exacer + race + unemployed + age.grp
    #           + calg1 + cs14 + cs16 + epsmean + qoltot + Chg.pansstotal + pct.gain + phase.change.cum.rec
    #           + calg1*calg1 + cs14*cs14 + cs16*cs16 + epsmean*epsmean + qoltot*qoltot + Chg.pansstotal*Chg.pansstotal + pct.gain*pct.gain + phase.change.cum.rec
    #           , family=binomial("logit")
    #)
    
    fit   <- gam(data=indata.fit,
                 studydisc ~
                 + s(time,4) 
                 + s(Bpansstotal,4) 
                 + site.ro+site.sh+site.uc+site.va
                 + exacer + race + unemployed + age.grp
                 + s(calg1,4) + s(cs14,4) + s(cs16,4) + s(epsmean,4) + s(qoltot,4) + s(Chg.pansstotal,4) + s(pct.gain,4) + phase.change.cum.rec
                 + s(cs14*time,4) + s(cs14*time,4) + s(cs16*time,4) + s(epsmean*time,4) + s(qoltot*time,4) + s(Chg.pansstotal*time,4) + s(pct.gain*time,4) + s(phase.change.cum.rec*time,4)  
                 ,family=binomial("logit")
                 )

    print(fit)
    
    indata.pred <- indata %>% filter(source=="pred")
    
    X1 <- predict(fit,newdata=indata.pred,type="response")
    X0 <- 1-X1
    pmat <- data.frame(X1,X0)
    names(pmat) <- c("1","0")
    output <- rMultinom(pmat,1) %>% data.frame
    
  }
  
  result <- fxn.data %>% group_by(!!! group_var) %>% do(mod.pred(.)) %>% ungroup() 
  
}

temp <- pred.studydisc(catie,step6.trial,treat.grp)
colnames(temp) <- c("treat.grp","studydisc")

step7.trial <- data.frame(step6.trial,temp$studydisc) %>% rename(studydisc=temp.studydisc)

#######################
# STEP 8: format data #
#######################

step8.trial <- step7.trial %>%
  group_by(CATIEID) %>%
  arrange(CATIEID,time) %>%
  mutate(studydisc=ifelse(.data$studydisc=="1",1,
                                        ifelse(.data$studydisc=="0",0,NA)),
         studydisc=cumsum(.data$studydisc),
         studydisc=cumsum(.data$studydisc)) %>%
  filter(.data$studydisc<=1)

#######################
# check distributions #
#######################

catie %>% group_by(time) %>% select(pansstotal,cs14,cs16,weight,epsmean,qoltot,Chg.pansstotal,pct.gain,phase.change.cum.rec,studydisc) %>% summarize_all("mean",na.rm=TRUE)
step8.trial %>% group_by(time) %>% select(pansstotal,cs14,cs16,weight,epsmean,qoltot,Chg.pansstotal,pct.gain,phase.change.cum.rec,studydisc) %>% summarize_all("mean",na.rm=TRUE)

###################################################################
## STEP 9: estimate probabilities for inverse probability weights #
###################################################################

p.est <- function (input,strategy,response="studydisc") {
  
  #input     <- step8.trial
  #response  = "studydisc"
  #strategy  = "pool"
  
  #print model output
  print.model <- function (input) {
    
    coef <- input %>% coef() %>% as.vector() 
    se   <- input %>% confint.default() %>% as.matrix() 
    output <- cbind(coef,se)
    output <- exp(output) %>% round(2)
    output[,1:3] <- round(output[,1:3],2)  
    return(output)
  }
  
  indata <- input %>% rename(Y=response)
  
  mergedata <- indata
  indata    <- indata %>% filter(.data$time!=18) 
  outdata   <- indata
  
  ## BUILD THE MODELS ##
  
  #exposure
  indata.0 <- indata %>% ungroup() %>% filter(time==0)  
  
  num.x   <- vgam(formula=treat.grp~td+zprcort+td*zprcort
                  + s_(Bpansstotal,4)
                  +exacer + age.grp + race + educ.bin + unemployed
                  + site.ro+site.sh+site.uc+site.va
                  ,family=multinomial(),data=indata.0)
  
  den.x.b <- vgam(formula=treat.grp~td+zprcort+td*zprcort
                  +exacer + s_(Bpansstotal,4) + age.grp + race + educ.bin + unemployed
                  + site.ro+site.sh+site.uc+site.va
                  + s_(Bcalg1,4) + s_(Bcs14,4) + s_(Bcs16,4) + s_(Bepsmean,4) +  s_(Bqoltot,4)
                  ,family=multinomial(),data=indata.0)
  
  print("Treatment Numerator Model")
  print(num.x)
  print(print.model(num.x))
  print("Treatment Denominator Model")
  print(den.x.b)
  print(print.model(den.x.b))
  
  #dropout    
  if (strategy=="pool") {
    
    num.p <- gam(indata,
                 
                 formula = Y~
                   treat.grp
                 + s(time,4)  
                 + s(Bpansstotal,4)
                 + exacer + age.grp + race + educ.bin + unemployed
                 + site.ro+site.sh+site.uc+site.va#+site.mc+site.pn+site.pp
                 ,family=binomial("logit") )   
    
    den.p.b  <- gam(indata,
                    formula = Y ~ 
                      treat.grp
                    + s(time,4) 
                    + s(Bpansstotal,4) 
                    + site.ro+site.sh+site.uc+site.va#+site.mc+site.pn+site.pp
                    + exacer + race + unemployed + age.grp
                    + s(Bcalg1) + s(Bcs14,4) + s(Bcs16,4) + s(Bepsmean,4) + s(Bqoltot,4)  
                    + s(Bcalg1*time,4) + s(Bcs14*time,4) + s(Bcs16*time,4) + s(Bepsmean*time,4) + s(Bqoltot*time,4)  
                    ,family=binomial("logit"))  
    
    den.p.t   <- gam(indata,
                     formula = Y ~ 
                       treat.grp
                     + s(time,4) 
                     + s(Bpansstotal,4) 
                     +site.ro+site.sh+site.uc+site.va#+site.mc+site.pn+site.pp
                     + exacer + race + unemployed + age.grp
                     + s(calg1,4) + s(cs14,4) + s(cs16,4) + s(epsmean,4) + s(qoltot,4) + s(Chg.pansstotal,4) + s(pct.gain,4) + phase.change.cum.rec
                     + s(cs14*time,4) + s(cs14*time,4) + s(cs16*time,4) + s(epsmean*time,4) + s(qoltot*time,4) + s(Chg.pansstotal*time,4) + s(pct.gain*time,4) + s(phase.change.cum.rec*time,4)  
                     ,family=binomial("logit"))
    
    print("Pooled Dropout Numerator Model")
    print(num.p)
    print(print.model(num.p))
    print("Pooled Dropout Denominator Model: Baseline Covs")
    print(den.p.b)
    print(print.model(den.p.b))
    print("Pooled Dropout Denominator Model: Time-varying Covs")
    print(den.p.t)
    print(print.model(den.p.t))
    
  } else if (strategy=="treat.specific") {
    
    fun.treat.num <- function (indata) {
      outdata <- indata 
      fit <- gam(indata,
                 formula = Y ~ 
                   + s(time,4) 
                 + s(Bpansstotal,4)
                 +exacer + age.grp + race + educ.bin + unemployed
                 +site.ro+site.sh+site.uc+site.va#+site.mc+site.pn+site.pp
                 ,family=binomial("logit") )
      print("Treatment-specific Numerator model i: 1=OLA, 2=QUE, 3=RIS, 4=PER, 5=ZIP")
      print(fit)
      print(print.model(fit))
      outdata$pred <- predict(fit,indata,type="response") 
      outdata %>% ungroup() %>% select(CATIEID,time,pred) %>% rename(num.p=pred) %>% data.frame()
    }
    treat.pred.num <- indata %>% group_by(treat.grp) %>% do(fun.treat.num(.))
    
    fun.treat.den.b <- function (indata) {
      outdata <- indata 
      fit <- gam(indata,
                 formula = Y ~  
                   + s(time,4) 
                 + s(Bpansstotal,4) 
                 + site.ro+site.sh+site.uc+site.va#+site.mc+site.pn+site.pp
                 + exacer + race + unemployed + age.grp
                 +s(Bcalg1,4) + s(Bcs14,4) + s(Bcs16,4) + s(Bepsmean,4) + s(Bqoltot,4)  
                 + s(Bcalg1,4)*s(time,4) + s(Bcs14*time,4) + s(Bcs16*time,4) + s(Bepsmean*time,4) + s(Bqoltot*time,4)  
                 ,family=binomial("logit"))
      print("Treatment-specific Denominator, Baseline Covs model i: 1=OLA, 2=QUE, 3=RIS, 4=PER, 5=ZIP")
      print(fit)
      print(print.model(fit))
      outdata$pred <- predict(fit,indata,type="response") 
      outdata %>% ungroup() %>% select(CATIEID,time,pred) %>% rename(den.p.b=pred) %>% data.frame()
    }
    treat.pred.den.b <- indata %>% group_by(treat.grp) %>% do(fun.treat.den.b(.))
    
    fun.treat.den.t <- function (indata) {
      outdata <- indata 
      fit <- gam(indata,
                 formula = Y ~ 
                   + s(time,4) 
                 + s(Bpansstotal,4) 
                 +site.ro+site.sh+site.uc+site.va#+site.mc+site.pn+site.pp
                 + exacer + race + unemployed + age.grp
                 + s(calg1,4) + s(cs14,4) + s(cs16,4) + s(epsmean,4) + s(qoltot,4)  + s(Chg.pansstotal,4) + s(pct.gain,4) + phase.change.cum.rec  
                 + s(cs14*time,4) + s(cs16*time,4) + s(epsmean*time,4) + s(qoltot*time,4) + s(Chg.pansstotal*time,4) + s(pct.gain*time,4) + s(phase.change.cum.rec*time,4)
                 ,family=binomial("logit"))
      print("Treatment-specific Denominator, Time-varying Covs model i: 1=OLA, 2=QUE, 3=RIS, 4=PER, 5=ZIP")
      print(fit)
      print(print.model(fit))
      outdata$pred <- predict(fit,indata,type="response") 
      outdata %>% ungroup() %>% select(CATIEID,time,pred) %>% rename(den.p.t=pred) %>% data.frame()
    }
    treat.pred.den.t <- indata %>% group_by(treat.grp) %>% do(fun.treat.den.t(.))
    
  }
  
  ## PREDICT PROBABILITIES ##
  
  indata[,c(c("num.ola","num.que","num.ris","num.per","num.zip"),c("den.ola","den.que","den.ris","den.per","den.zip"))] <- NA
  outdata[,c(c("num.ola","num.que","num.ris","num.per","num.zip"),c("den.ola","den.que","den.ris","den.per","den.zip"))] <- NA
  
  outdata[which(outdata$time==0),c("num.ola","num.que","num.ris","num.per","num.zip")] <- predict(object=num.x  ,newdata=indata.0,type="response")
  outdata[which(outdata$time==0),c("den.ola","den.que","den.ris","den.per","den.zip")] <- predict(object=den.x.b,newdata=indata.0,type="response")
  
  if (strategy=="pool") {
    
    outdata[,"num.p"]   <- predict(object=num.p,newdata=indata,type="response")
    outdata[,"den.p.b"] <- predict(object=den.p.b,newdata=indata,type="response")
    outdata[,"den.p.t"] <- predict(object=den.p.t,newdata=indata,type="response")
    
  } else if (strategy=="treat.specific") {
    
    outdata.a     <- outdata %>% ungroup() %>% arrange(CATIEID,time)
    outdata.b     <- left_join(outdata.a,treat.pred.num,by=c("CATIEID","time"))
    outdata.c     <- left_join(outdata.b,treat.pred.den.b,by=c("CATIEID","time"))
    outdata       <- left_join(outdata.c,treat.pred.den.t,by=c("CATIEID","time")) 
    
  }
  
  outdata <- outdata %>% 
    arrange(CATIEID,time) %>% 
    select(CATIEID,time,num.p,den.p.b,den.p.t,num.ola,num.que,num.ris,num.per,num.zip,den.ola,den.que,den.ris,den.per,den.zip,Y) 
  
  mergedata <- mergedata %>%
    arrange(CATIEID,time) 
  
  finaldata <- left_join(mergedata,outdata,by=c("CATIEID","time","Y")) 
  finaldata[which(finaldata$time==18 & finaldata$Y==0),c("num.p","den.p.b","den.p.t","num.ola","num.que","num.ris","num.per","num.zip","den.ola","den.que","den.ris","den.per","den.zip")] <- 0
  finaldata[which(finaldata$time==18 & finaldata$Y==1),c("num.p","den.p.b","den.p.t","num.ola","num.que","num.ris","num.per","num.zip","den.ola","den.que","den.ris","den.per","den.zip")] <- 1
  
  if (strategy=="pool") {
    
    finaldata <- rename(finaldata,num.po=num.p,den.po.b=den.p.b,den.po.t=den.p.t)
    
  } else if (strategy=="treat.specific") {
    
    finaldata <- rename(finaldata,num.tr=num.p,den.tr.b=den.p.b,den.tr.t=den.p.t)
    
  } 
  
  s_response <- sym(response)
  
  finaldata <- finaldata %>%
    rename(!! s_response := Y) %>%
    data.frame()  
  
}

#estimated predicted probabilities of treatment assignment and dropout

pred.pool  <- step8.trial %>% p.est(strategy="pool")	
pred.treat <- step8.trial %>% p.est(strategy="treat.specific")

#merge in predicted probabilities of treatment assignment and dropout
step9.trial.pool  <- pred.pool  %>% select(c(CATIEID,time,starts_with("num"),starts_with("den")))   
step9.trial.treat <- pred.treat %>% select(c(CATIEID,time,num.tr,den.tr.b,den.tr.t))   
temp              <- step8.trial %>% ungroup() %>% arrange(CATIEID,time)
temp.po           <- left_join(temp,step9.trial.pool,by=c("CATIEID","time"))
temp.tr           <- left_join(temp.po,step9.trial.treat,by=c("CATIEID","time"))
step9.trial       <- temp.tr

#################################################
## STEP 10: create inverse probability weights ##
#################################################

#FUNCTION TO CREATE UNSTABLISZED & STABILIZED TIME-SPECIFIC, CUMULATIVE WEIGHTS FOR STUDY DROPOUT
makeweight <- function(input,truncate,limit) {
  
  interim <- input %>% arrange(CATIEID,time) %>% mutate(
    
    #dropout indicator
    s=as.numeric(studydisc),    
    
    #probabilities of treatment
    num.x=ifelse(time==0,(treat.grp.ola*(num.ola)
                          +treat.grp.que*(num.que)
                          +treat.grp.ris*(num.ris)
                          +treat.grp.per*(num.per)
                          +treat.grp.zip*(num.zip)),1),
    
    num.x=cumprod(num.x),

    den.x=ifelse(time==0,(treat.grp.ola*(den.ola)
                            +treat.grp.que*(den.que)
                            +treat.grp.ris*(den.ris)
                            +treat.grp.per*(den.per)
                            +treat.grp.zip*(den.zip)),1),
    
    den.x=cumprod(den.x),
    
    
    #unstabilized treatment weights
    uwx.b=ifelse(time==0,1/(treat.grp.ola*(den.ola)
                            +treat.grp.que*(den.que)
                            +treat.grp.ris*(den.ris)
                            +treat.grp.per*(den.per)
                            +treat.grp.zip*(den.zip)),1),
    
    #unstabilized dropout weights
    uwpo.b=1/(s*den.po.b+(1-s)*(1-den.po.b)),
    uwpo.t=1/(s*den.po.t+(1-s)*(1-den.po.t)),
    uwtr.b=1/(s*den.tr.b+(1-s)*(1-den.tr.b)),
    uwtr.t=1/(s*den.tr.t+(1-s)*(1-den.tr.t)),
    
    #stabilized treatment weights
    wx.b=ifelse(time==0,(treat.grp.ola*(num.ola)
                         +treat.grp.que*(num.que)
                         +treat.grp.ris*(num.ris)
                         +treat.grp.per*(num.per)
                         +treat.grp.zip*(num.zip))/(treat.grp.ola*(den.ola)
                                                    +treat.grp.que*(den.que)
                                                    +treat.grp.ris*(den.ris)
                                                    +treat.grp.per*(den.per)
                                                    +treat.grp.zip*(den.zip)),1),
    #stabilized dropout 
    wpo.b=(s*num.po+(1-s)*(1-num.po))/(s*den.po.b+(1-s)*(1-den.po.b)),
    wpo.t=(s*num.po+(1-s)*(1-num.po))/(s*den.po.t+(1-s)*(1-den.po.t)),
    wtr.b=(s*num.tr+(1-s)*(1-num.tr))/(s*den.tr.b+(1-s)*(1-den.tr.b)),
    wtr.t=(s*num.tr+(1-s)*(1-num.tr))/(s*den.tr.t+(1-s)*(1-den.tr.t))
  )
  
  if(truncate=="yes") {
    
    interim <- interim %>% mutate(
      
      uwx.b=ifelse(uwx.b>quantile(uwx.b,limit,na.rm=TRUE),quantile(uwx.b,limit,na.rm=TRUE),uwx.b),
      uwpo.b=ifelse(uwpo.b>quantile(uwpo.b,limit,na.rm=TRUE),quantile(uwpo.b,limit,na.rm=TRUE),uwpo.b),
      uwpo.t=ifelse(uwpo.t>quantile(uwpo.t,limit,na.rm=TRUE),quantile(uwpo.t,limit,na.rm=TRUE),uwpo.t),
      uwtr.b=ifelse(uwtr.b>quantile(uwtr.b,limit,na.rm=TRUE),quantile(uwtr.b,limit,na.rm=TRUE),uwtr.b),
      uwtr.t=ifelse(uwtr.t>quantile(uwtr.t,limit,na.rm=TRUE),quantile(uwtr.t,limit,na.rm=TRUE),uwtr.t),
      
      wx.b=ifelse(wx.b>quantile(wx.b,limit,na.rm=TRUE),quantile(wx.b,limit,na.rm=TRUE),wx.b),
      wpo.b=ifelse(wpo.b>quantile(wpo.b,limit,na.rm=TRUE),quantile(wpo.b,limit,na.rm=TRUE),wpo.b),
      wpo.t=ifelse(wpo.t>quantile(wpo.t,limit,na.rm=TRUE),quantile(wpo.t,limit,na.rm=TRUE),wpo.t),
      wtr.b=ifelse(wtr.b>quantile(wtr.b,limit,na.rm=TRUE),quantile(wtr.b,limit,na.rm=TRUE),wtr.b),
      wtr.t=ifelse(wtr.t>quantile(wtr.t,limit,na.rm=TRUE),quantile(wtr.t,limit,na.rm=TRUE),wtr.t)
      
    )
    
  }
  
  ## CALCULATE CUMULATIVE WEIGHTS ##  
  
  output <- interim %>% 
    arrange(CATIEID,time) %>% 
    group_by(CATIEID) %>% 
    mutate(cum.uwx.b=cumprod(uwx.b),
           cum.uwpo.b=cum.uwx.b*cumprod(uwpo.b),
           cum.uwpo.t=cum.uwx.b*cumprod(uwpo.t),
           cum.uwtr.b=cum.uwx.b*cumprod(uwtr.b),
           cum.uwtr.t=cum.uwx.b*cumprod(uwtr.t),
           
           cum.wx.b=cumprod(wx.b),
           cum.wpo.b=cum.wx.b*cumprod(wpo.b),
           cum.wpo.t=cum.wx.b*cumprod(wpo.t),
           cum.wtr.b=cum.wx.b*cumprod(wtr.b),
           cum.wtr.t=cum.wx.b*cumprod(wtr.t)
    )
  
}								  

step9.trial.simNo <- makeweight(step9.trial,truncate="no") %>% select(CATIEID,time,num.x,den.x,wx.b,num.po,den.po.t,num.tr,den.tr.t,cum.wpo.t,cum.wtr.t) %>% rename(wpo=cum.wpo.t,wtr=cum.wtr.t)
step9.trial.sim99 <- makeweight(step9.trial,truncate="yes",limit=.99) %>% select(CATIEID,time,cum.wpo.t,cum.wtr.t) %>% rename(wpo99=cum.wpo.t,wtr99=cum.wtr.t) 
step9.trial.sim95 <- makeweight(step9.trial,truncate="yes",limit=.95) %>% select(CATIEID,time,cum.wpo.t,cum.wtr.t) %>% rename(wpo95=cum.wpo.t,wtr95=cum.wtr.t)
step9.trial.sim90 <- makeweight(step9.trial,truncate="yes",limit=.90) %>% select(CATIEID,time,cum.wpo.t,cum.wtr.t) %>% rename(wpo90=cum.wpo.t,wtr90=cum.wtr.t)

step10.trial <- step9.trial %>%
  select(CATIEID:studydisc) %>%
  left_join(step9.trial.simNo,by=c("CATIEID","time")) %>%
  left_join(step9.trial.sim99,by=c("CATIEID","time")) %>%
  left_join(step9.trial.sim95,by=c("CATIEID","time")) %>%
  left_join(step9.trial.sim90,by=c("CATIEID","time"))

####################
## MASK ID NUMBERS #
####################

catie.sim <- step10.trial %>%
  mutate(CATIEID=(CATIEID-1990)*(-1)+1990-81)

#############
# SAVE DATA #
#############

write_csv(catie.sim,paste(path.data,"catie_sim.csv",sep=""))
