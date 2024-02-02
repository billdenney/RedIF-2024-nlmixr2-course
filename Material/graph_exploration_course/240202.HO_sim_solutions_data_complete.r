
# Libraries ---------------------------------------------------------------

rm(list=ls())


library(tidyr)
library(SimTools)
library(dplyr)

library(survminer)
library(survival)

library(MASS)
library(Hmisc)
library(gam)

library(PKNCA)


# PKPD simulations settings brainstorm -----------------------------------------------------

# fictive PD

# Safety driven by Cplasma -> example albumin drop
log(2)/0.0015/24 # 19d half life: kout = 0.0015 h-1 in human
40 * 0.0015 # 40 g/L normal levels: kin = 0.06 g/L/h
# Our effect will decrease Kin -> marker of liver damage

# species differences
# Note: SPEC controls typical body weight; 1=human,2=monkey,3=rat,4=mouse
# For mab, add nonlinear cl (reflective of some "non brain" TMDD), for the MONKEY but not RAT (since rat less X react)
# For mab, add ADA issues, for rat: 50%, monkey 10%

# PKPD - Simulate non clinical ------------------------------

# User input
source("Translation_PKPD.r")
nid <- 10
prob_monkey_ada <- 0.1
prob_rodent_ada <- 0.5
ada_onset <- 3*24 

# Fictive Safety and Efficacy Turnover
# Note: I call kout TV, because in the model I'll scale with body weight (-0.25), not for Kin (choice at this point)
# Note: I call kin TV, because in the model we'll add IIV
TV_SAF_KOUT<-0.0015 ; TV_SAF_KIN<-40*TV_SAF_KOUT
TV_EFF_KOUT<-0.15 ; TV_EFF_KIN<-10*TV_EFF_KOUT

# Fictive mAb (based on Betts et al); CL mL/h/kg 0.15 -> 0.15*70/1000;V1 mL/kg 46.31;Q mL/h/kg 0.27;V2 mL/kg 31.47
p_mab <- c(TVCL=0.01,TVV1=3,TVKA=0.1,TVF=0.8,TVQ=0.02,TVV2=3,TVVMAX=1,TVKM=0.1,
           KP = 0.01 , TVQB=0.02, TVVB = 0.03,
           SAF_MAX =    0, SAF_C50 = 10, SAF_GAM = 2, EFF_MAX = -1, EFF_C50 = 2, EFF_GAM = 1) #Volume L; Time h; Dose mg (conc mg/L)
# Fictive small molecule
p_sm  <- c(TVCL=0.4,TVV1=3,TVKA=0.05,TVF=0.6,TVQ=0.8,TVV2=9,TVVMAX=0,TVKM=1,
           KP = 1    , TVQB=0.02, TVVB = 0.03,
           SAF_MAX = -0.5, SAF_C50 = 10, SAF_GAM = 2, EFF_MAX = -1, EFF_C50 = 2, EFF_GAM = 1) #Volume L; Time h

# IIV
var1 <- 0.1     # variance CL
var2 <- 0.04    # variance V1
var3 <- 0.1     # variance VMAX
var4 <- 0.1     # variance KA
var5 <- 0.1     # variance KIN SAFETY
var6 <- 0.1     # variance KIN EFFICACY

# Error
sd_cplasma_prp_err <- 0.3
sd_cbrain_prp_err <- 0.4
sd_eff_prp_err <- 0.15
sd_saf_prp_err <- 0.1

# LLOQ
lloq_cplasma <- 0.05
lloq_cbrain <- 0.05
lloq_eff <- 0.05
lloq_saf <- 0.05

s_df <- expand.grid(SPEC = c(2,3,4),
                    DOSE_MGKG = c(3,10,30),
                    ROUTE = c(1,2),
                    II = c(48),
                    NDOS = c(1,10),
                    CMPD = c(1,2)
                    ) %>%
  mutate(
         TVCL=ifelse(CMPD==1,p_mab[["TVCL"]],p_sm[["TVCL"]]),
         TVV1=ifelse(CMPD==1,p_mab[["TVV1"]],p_sm[["TVV1"]]),
         TVKA=ifelse(CMPD==1,p_mab[["TVKA"]],p_sm[["TVKA"]]),
         TVF=ifelse(CMPD==1,p_mab[["TVF"]],p_sm[["TVF"]]),
         TVQ=ifelse(CMPD==1,p_mab[["TVQ"]],p_sm[["TVQ"]]),
         TVV2=ifelse(CMPD==1,p_mab[["TVV2"]],p_sm[["TVV2"]]),
         TVVMAX=ifelse(CMPD==1,p_mab[["TVVMAX"]],p_sm[["TVVMAX"]]),
         TVVMAX=ifelse(CMPD==1&SPEC%in%c(3,4),0*TVVMAX,TVVMAX),        # Rodent mab non linearity not picked up
         TVKM=ifelse(CMPD==1,p_mab[["TVKM"]],p_sm[["TVKM"]]),
         KP=ifelse(CMPD==1,p_mab[["KP"]],p_sm[["KP"]]),
         TVQB=ifelse(CMPD==1,p_mab[["TVQB"]],p_sm[["TVQB"]]),
         TVVB=ifelse(CMPD==1,p_mab[["TVVB"]],p_sm[["TVVB"]]),
         TV_SAF_KOUT=TV_SAF_KOUT,TV_SAF_KIN=TV_SAF_KIN,TV_EFF_KOUT=TV_EFF_KOUT,TV_EFF_KIN=TV_EFF_KIN,
         SAF_MAX=ifelse(CMPD==1,p_mab[["SAF_MAX"]],p_sm[["SAF_MAX"]]),
         SAF_C50=ifelse(CMPD==1,p_mab[["SAF_C50"]],p_sm[["SAF_C50"]]),
         SAF_GAM=ifelse(CMPD==1,p_mab[["SAF_GAM"]],p_sm[["SAF_GAM"]]),
         EFF_MAX=ifelse(CMPD==1,p_mab[["EFF_MAX"]],p_sm[["EFF_MAX"]]),
         EFF_C50=ifelse(CMPD==1,p_mab[["EFF_C50"]],p_sm[["EFF_C50"]]),
         EFF_C50=ifelse(CMPD==1&SPEC%in%c(3,4),10*EFF_C50,EFF_C50),   # Rodent EC50 higher for mab
         EFF_GAM=ifelse(CMPD==1,p_mab[["EFF_GAM"]],p_sm[["EFF_GAM"]]),
         BW=ifelse(SPEC==1,70,ifelse(SPEC==2,3,ifelse(SPEC==3,0.3,ifelse(SPEC==4,0.02,-999)))), # -999 to ensure error as not wanted input
         DOSE_MG=DOSE_MGKG*BW)
s_df$SCEN <- 1:nrow(s_df)
nrep        <- nrow(s_df)

# eta generation
samp_shell <- data.frame(ID= 1:(nid*nrep), SCEN = rep(1:nrep,nid), REP = rep(1:nid,each=nrep))
set.seed(1) ; eta1 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var1))
set.seed(2) ; eta2 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var2))
set.seed(3) ; eta3 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var3))
set.seed(4) ; eta4 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var4))
set.seed(5) ; eta5 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var5))
set.seed(6) ; eta6 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var6))
samp_ETA <-  mutate(samp_shell,ETA1=eta1,ETA2=eta2,ETA3=eta3,ETA4=eta4,ETA5=eta5,ETA6=eta6)

set.seed(10) ; adadraw <- runif(nrow(samp_ETA))
samp <- dplyr::left_join(samp_ETA,s_df,by="SCEN") %>%
  mutate(ADADRAW=adadraw,ADA=ifelse(CMPD==1&SPEC==2&ADADRAW<prob_monkey_ada,1,ifelse(CMPD==1&SPEC%in%c(3,4)&ADADRAW<prob_rodent_ada,1,0)))

# Testing via # x <- 1 # unlist(samp[x,])
out <- lapply(1:nrow(samp), function(x){
  print(x)
  parm  <- par_func(TVCL=samp[x,]$TVCL,TVV1=samp[x,]$TVV1,TVQ=samp[x,]$TVQ,TVV2=samp[x,]$TVV2,TVVMAX=samp[x,]$TVVMAX,TVKM=samp[x,]$TVKM,TVKA=samp[x,]$TVKA,
                    ETA1=samp[x,]$ETA1,ETA2=samp[x,]$ETA2,ETA3=samp[x,]$ETA3,ETA4=samp[x,]$ETA4,ETA5=samp[x,]$ETA5,ETA6=samp[x,]$ETA6,
                    KP=samp[x,]$KP,TVQB=samp[x,]$TVQB,TVVB=samp[x,]$TVVB,
                    TV_SAF_KOUT=samp[x,]$TV_SAF_KOUT,TV_SAF_KIN=samp[x,]$TV_SAF_KIN,TV_EFF_KOUT=samp[x,]$TV_EFF_KOUT,TV_EFF_KIN=samp[x,]$TV_EFF_KIN,
                    SAF_MAX=samp[x,]$SAF_MAX,SAF_C50=samp[x,]$SAF_C50,SAF_GAM=samp[x,]$SAF_GAM,
                    EFF_MAX=samp[x,]$EFF_MAX,EFF_C50=samp[x,]$EFF_C50,EFF_GAM=samp[x,]$EFF_GAM,
                    SPEC=samp[x,]$SPEC,
                    BW=samp[x,]$BW,
                    ADA=samp[x,]$ADA,
                    ADASTART=ada_onset) # If ADA=1, ADA effect kicks in onset of time, currently no IIV
  inits <- init_func(parm)
  if(unlist(samp[x,])[["ROUTE"]]==1) evnt  <- dose_func(cmt=1,
                                                        value=samp[x,]$DOSE_MG,
                                                        # tinf=2, # careful for infusion dose in cmt=3, for bolus dose in cmt=1
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  if(unlist(samp[x,])[["ROUTE"]]!=1) evnt  <- dose_func(cmt=5,
                                                        value=samp[x,]$DOSE_MG * samp[x,]$TVF,
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  
  times  <- sort(c(0.083,0.25,0.5,1,2,seq(0,2*7*24,4)))
  data.frame(deSolve::lsoda(inits,times,des_func,parm,events=list(data=evnt)),
             ID=unlist(samp[x,])[["ID"]]
  )
})
out <- do.call(rbind,out)

outsel <- dplyr::select(out,time,C1,AUC=A4,CB,SAFETY=A7,EFFICACY=A8,ID) %>% 
  dplyr::left_join(dplyr::select(samp,ID:REP,SPEC:CMPD,BW,DOSE_MG,ADA),by="ID") %>%
  mutate(SPEC=ifelse(SPEC==1,"Human",ifelse(SPEC==2,"Monkey",ifelse(SPEC==3,"Rat",ifelse(SPEC==4,"Mouse","ERROR")))),
         CMPD=as.factor(ifelse(CMPD==1,"mAb","SM")))

outsel0 <- subset(outsel,time==0) %>%
  mutate(SAFETY0 = SAFETY, EFFICACY0 = EFFICACY) %>%
  dplyr::select(ID,SAFETY0,EFFICACY0)

outsel <- dplyr::left_join(outsel,outsel0)

# Add error
set.seed(10) ; p_err_cplasma <- rnorm(n = nrow(outsel), mean = 0, sd = sd_cplasma_prp_err)
set.seed(20) ; p_err_cbrain <- rnorm(n = nrow(outsel), mean = 0, sd = sd_cbrain_prp_err)
set.seed(30) ; p_err_eff <- rnorm(n = nrow(outsel), mean = 0, sd = sd_eff_prp_err)
set.seed(40) ; p_err_saf <- rnorm(n = nrow(outsel), mean = 0, sd = sd_saf_prp_err)

outsel_err <- mutate(outsel,
                     C1_ERR = p_err_cplasma, C1_Y_ORI = C1*(1+C1_ERR), C1_BQL = ifelse(C1_Y_ORI<lloq_cplasma,1,0), C1_Y = ifelse(C1_BQL==1,0,C1_Y_ORI),
                     CB_ERR = p_err_cbrain, CB_Y_ORI = CB*(1+CB_ERR), CB_BQL = ifelse(CB_Y_ORI<lloq_cbrain,1,0), CB_Y = ifelse(CB_BQL==1,0,CB_Y_ORI),
                     EFF_ERR = p_err_eff, EFF_Y_ORI = EFFICACY*(1+EFF_ERR), EFF_BQL = ifelse(EFF_Y_ORI<lloq_eff,1,0), EFF_Y = ifelse(EFF_BQL==1,0,EFF_Y_ORI),
                     SAF_ERR = p_err_saf, SAF_Y_ORI = SAFETY*(1+SAF_ERR), SAF_BQL = ifelse(SAF_Y_ORI<lloq_saf,1,0), SAF_Y = ifelse(SAF_BQL==1,0,SAF_Y_ORI)) 

write.csv(outsel_err, file=paste0("outsel_err","_n",nid,".csv"), row.names = FALSE, quote = FALSE, na = ".")

# Ipred Explore Step 1-2-3-4 ----------------------------------------------

# Explore Step 1 SM IV PK sd - IPRED

ggplot(subset(outsel,time<7*24&ROUTE==1&NDOS==1&CMPD=="SM"&SPEC%in%c("Rat")))+
  geom_line(aes(x=time,y=C1,group=ID,col=as.factor(ROUTE)))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_y_log10()+
  geom_hline(yintercept = 0.05)+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

ggplot(subset(outsel,time<7*24&ROUTE==1&NDOS==1&CMPD=="SM"))+
  geom_line(aes(x=time,y=C1,group=ID,col=as.factor(ROUTE)))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_y_log10()+
  geom_hline(yintercept = 0.05)+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

# Explore Step 2 SM IV/PO PK sd - IPRED

ggplot(subset(outsel,time<7*24&NDOS==1&CMPD=="SM"&SPEC%in%c("Rat")))+
  geom_line(aes(x=time,y=C1,group=ID,col=as.factor(ROUTE)))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_y_log10()+
  geom_hline(yintercept = 0.05)+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

# Explore Step 3 mAb PK sd - IPRED

ggplot(subset(outsel,time<14*24&NDOS==1&CMPD=="mAb"&SPEC%in%c("Monkey","Rat")&ROUTE==1))+
  geom_line(aes(x=time,y=C1,group=ID,col=CMPD))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_y_log10()+
  geom_hline(yintercept = c(0.05,1))+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

# Explore Step 4 PKPD q2d - IPRED
ggplot(subset(outsel,SPEC%in%c("Monkey","Rat")&NDOS>1&ROUTE==2))+
  geom_line(aes(x=time,y=C1,group=ID,col=CMPD))+
  scale_y_log10()+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  geom_hline(yintercept = 10)+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

ggplot(subset(outsel,SPEC%in%c("Monkey","Rat")&NDOS>1&ROUTE==2))+
  geom_line(aes(x=time,y=SAFETY,group=ID,col=CMPD))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

ggplot(subset(outsel,SPEC%in%c("Monkey","Rat")&NDOS>1&ROUTE==2))+
  geom_line(aes(x=time,y=100*(SAFETY-SAFETY0)/SAFETY0,group=ID,col=CMPD))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

ggplot(subset(outsel,SPEC%in%c("Monkey","Rat")&NDOS>1&ROUTE==2))+
  geom_line(aes(x=time,y=CB,group=ID,col=CMPD))+
  scale_y_log10()+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  geom_hline(yintercept = 2)+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

ggplot(subset(outsel,SPEC%in%c("Monkey","Rat")&NDOS>1&ROUTE==2))+
  geom_line(aes(x=time,y=EFFICACY,group=ID,col=CMPD))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

ggplot(subset(outsel,SPEC%in%c("Monkey","Rat")&NDOS>1&ROUTE==2))+
  geom_line(aes(x=time,y=100*(EFFICACY-EFFICACY0)/EFFICACY0,group=ID,col=CMPD))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

# Step 1-2-3: SD PK - SELECT TPs ------------------------------
trans_sdc1_err <- read.csv("outsel_err_n10.csv") %>%
  subset(NDOS==1&CMPD=="mAb"&time%in%c(0.083,0.5,1,2,4,24,24*seq(2,14,4)))

trans_sdc2_err <- read.csv("outsel_err_n10.csv") %>%
  subset(NDOS==1&CMPD=="SM"&time%in%c(0.083,0.5,1,2,4,24,24*seq(3,7,2)))

trans_sd_err <- rbind(trans_sdc1_err,trans_sdc2_err)

write.csv(subset(dplyr::select(trans_sd_err,ID,TIME=time,SPEC,DOSE_MGKG,ROUTE,CMPD,BW,ADA,Cplasma=C1_Y,BQL=C1_BQL),
                 CMPD=="SM"&SPEC=="Rat"&ROUTE==1), 
          file="./HO/Step1/sm_rat_pk_iv_sd.csv", row.names = FALSE, quote = FALSE)

write.csv(subset(dplyr::select(trans_sd_err,ID,TIME=time,SPEC,DOSE_MGKG,ROUTE,CMPD,BW,ADA,Cplasma=C1_Y,BQL=C1_BQL),
                 CMPD=="SM"&ROUTE==1), 
          file="./HO/Step1/sm_trans_pk_iv_sd.csv", row.names = FALSE, quote = FALSE)

write.csv(subset(dplyr::select(trans_sd_err,ID,TIME=time,SPEC,DOSE_MGKG,ROUTE,CMPD,BW,ADA,Cplasma=C1_Y,BQL=C1_BQL),
                 CMPD=="SM"&SPEC=="Rat"), 
          file="./HO/Step2/sm_rat_pk_ivpo_sd.csv", row.names = FALSE, quote = FALSE)

write.csv(subset(dplyr::select(trans_sd_err,ID,TIME=time,SPEC,DOSE_MGKG,ROUTE,CMPD,BW,ADA,Cplasma=C1_Y,BQL=C1_BQL),
                 CMPD=="mAb"&ROUTE==1&SPEC%in%c("Rat","Monkey")), 
          file="./HO/Step3/mab_trans_pk_iv_sd.csv", row.names = FALSE, quote = FALSE)

# Step 1-2-3: SD PK - NCA ------------------------------

d_conc <- subset(trans_sd_err) %>% dplyr::select(Time=time,conc=C1_Y,Subject=ID,Dose=DOSE_MG)
d_dose <- subset(d_conc,!duplicated(Subject)) %>% mutate(Time=0,conc=0)
my_start_time <- 0
d_conc <- rbind(d_conc,mutate(d_dose,Time=my_start_time,conc=0)) %>% arrange(Subject,Time)
head(d_conc)

conc_obj <-PKNCAconc(d_conc,conc~Time|Subject)
dose_obj <-PKNCAdose(d_dose,Dose~Time|Subject)

PKNCA.options() # look here what options to set on in intervals_manual, but could also change default here
# to reset 
# PKNCA.options(default=TRUE)

intervals_manual  <-
  data.frame(start=my_start_time, end=c(48, Inf),
             cmax=c(FALSE, TRUE),
             tmax=c(FALSE, TRUE),
             auclast=c(TRUE, FALSE),
             aucinf.obs=FALSE,
             aucinf.pred=c(FALSE, TRUE),
             cl.pred=c(FALSE, TRUE),
             vd.pred=c(FALSE, TRUE),
             vss.pred=c(FALSE, TRUE),
             aumcinf.pred=c(FALSE,TRUE))

data_obj <- PKNCAdata(conc_obj, dose_obj, intervals= intervals_manual)
results_obj <- pk.nca(data_obj)

nca_sd_trans_long <- results_obj$result %>% 
  subset(!PPTESTCD%in%c("span.ratio","r.squared")) %>%
  mutate(Parameter = ifelse(end=="Inf",PPTESTCD,paste0(PPTESTCD,"_0h_",end,"h"))) %>%
  dplyr::select(ID=Subject, Parameter ,Value=PPORRES) %>%
  dplyr::left_join(dplyr::select(subset(trans_sd_err,!duplicated(ID)),ID,CMPD,SPEC,BW,ROUTE,NDOS,DOSE_MGKG,ADA))

nca_sd_trans <- nca_sd_trans_long %>% tidyr::spread(key=Parameter, value=Value)
head(nca_sd_trans)

# Careful: depending on input of NCA, can be too sparse, will create NA's, which will make the input a character

write.csv(subset(nca_sd_trans,CMPD=="SM"&SPEC=="Rat"&ROUTE==1), 
          file="./HO/Step1/sm_rat_nca_iv_sd.csv", row.names = FALSE, quote = FALSE)

write.csv(subset(nca_sd_trans,CMPD=="SM"&ROUTE==1), 
          file="./HO/Step1/sm_trans_nca_iv_sd.csv", row.names = FALSE, quote = FALSE)

write.csv(subset(nca_sd_trans,CMPD=="SM"&SPEC=="Rat"), 
          file="./HO/Step2/sm_rat_nca_ivpo_sd.csv", row.names = FALSE, quote = FALSE)

write.csv(subset(nca_sd_trans,CMPD=="mAb"&ROUTE==1&SPEC%in%c("Rat","Monkey")), 
          file="./HO/Step3/mab_trans_nca_iv_sd.csv", row.names = FALSE, quote = FALSE)

# Step 1: Results - Learn allometry and dose linearity based on linear model ---------------------------------------------------------

pl_rat_pk_iv <- ggplot(subset(trans_sd_err,C1_BQL==0&ROUTE==1&CMPD=="SM"&SPEC%in%c("Rat")),
                       aes(x=time,y=C1_Y))+
  geom_point() +
  geom_line(aes(group=ID)) +
  # stat_summary(fun=median, geom='line', ) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,168,24))+
  facet_wrap(~DOSE_MGKG,labeller="label_both")+
  geom_hline(yintercept = 0.05)+
  labs(subtitle="iv sd SM rat",x="Time (h)",y="Cplasma (?g/mL)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[c(2,1)])+ theme(legend.position="bottom")
pl_rat_pk_iv


pl_rat_pk_iv_zoom <- ggplot(subset(trans_sd_err,C1_BQL==0&ROUTE==1&CMPD=="SM"&SPEC%in%c("Rat")),
                       aes(x=time,y=C1_Y))+
  geom_point() +
  geom_line(aes(group=ID)) +
  # stat_summary(fun=median, geom='line', ) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,168,2))+
  coord_cartesian(xlim=c(0,8))+
  facet_wrap(~DOSE_MGKG,labeller="label_both")+
  geom_hline(yintercept = 0.05)+
  labs(col="Route",subtitle="iv sd SM rat",x="Time (h)",y="Cplasma (?g/mL)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[c(2,1)])+ theme(legend.position="bottom")
pl_rat_pk_iv_zoom


allo_sm_trans_nca_iv_sd <- read.csv("./HO/Step1/sm_trans_nca_iv_sd.csv") %>%
  mutate(log.cl.pred=log(cl.pred),
         log.BW=log(BW),
         log.vd.pred=log(vd.pred),
         log.vss.pred=log(vss.pred),
         log.half.life=log(half.life))

pl_lin_auc <- ggplot(subset(allo_sm_trans_nca_iv_sd,SPEC=="Rat"),aes(x=DOSE_MGKG,y=aucinf.pred/DOSE_MGKG,group=as.factor(SPEC),col=as.factor(SPEC)))+
  geom_point() +
  scale_x_log10()+
  geom_smooth(method=lm , se=TRUE, formula = y ~ x) +
  labs(y="Dose normalized AUC (h*?g/mL / mg/kg)",x="Dose (mg/kg)",col="")+
  coord_cartesian(ylim=c(0,200))+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_lin_auc

mod.cl <- lm(log.cl.pred~log.BW,data=allo_sm_trans_nca_iv_sd)
mod.vd <- lm(log.vd.pred~log.BW,data=allo_sm_trans_nca_iv_sd)
mod.vss <- lm(log.vss.pred~log.BW,data=allo_sm_trans_nca_iv_sd)
mod.half.life <- lm(log.half.life~log.BW,data=allo_sm_trans_nca_iv_sd)

ci.cl <- paste0(signif(summary(mod.cl)$coef[[2,1]],2)," (",
                signif(summary(mod.cl)$coef[[2,1]]-1.96*summary(mod.cl)$coef[[2,2]],2),"-",
                signif(summary(mod.cl)$coef[[2,1]]+1.96*summary(mod.cl)$coef[[2,2]],2),")")

ci.vd <- paste0(signif(summary(mod.vd)$coef[[2,1]],2)," (",
                signif(summary(mod.vd)$coef[[2,1]]-1.96*summary(mod.vd)$coef[[2,2]],2),"-",
                signif(summary(mod.vd)$coef[[2,1]]+1.96*summary(mod.vd)$coef[[2,2]],2),")")

ci.vss <- paste0(signif(summary(mod.vss)$coef[[2,1]],2)," (",
                 signif(summary(mod.vss)$coef[[2,1]]-1.96*summary(mod.vss)$coef[[2,2]],2),"-",
                 signif(summary(mod.vss)$coef[[2,1]]+1.96*summary(mod.vss)$coef[[2,2]],2),")")

ci.half.life <- paste0(signif(summary(mod.half.life)$coef[[2,1]],2)," (",
                       signif(summary(mod.half.life)$coef[[2,1]]-1.96*summary(mod.half.life)$coef[[2,2]],2),"-",
                       signif(summary(mod.half.life)$coef[[2,1]]+1.96*summary(mod.half.life)$coef[[2,2]],2),")")

pl_lin_cl_notlog <- ggplot(allo_sm_trans_nca_iv_sd,aes(x=BW,y=cl.pred))+
  geom_point(aes(col=as.factor(DOSE_MGKG))) +
  labs(x="Body weight (kg)",col="Dose (mg/kg)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_lin_cl_notlog

pl_lin_cl <- ggplot(allo_sm_trans_nca_iv_sd,aes(x=BW,y=cl.pred))+
  geom_point(aes(col=as.factor(DOSE_MGKG))) +
  scale_x_log10()+scale_y_log10()+
  geom_smooth(method=lm , se=TRUE, formula = y ~ x,col="black") +
  labs(x="Body weight (kg)",col="Dose (mg/kg)",subtitle=paste0("Slope = Allometric coefficient = ",ci.cl))+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_lin_cl

ggplot(allo_sm_trans_nca_iv_sd,aes(x=vd.pred,y=vss.pred))+
  geom_point(aes(col=as.factor(DOSE_MGKG))) +
  scale_x_log10()+scale_y_log10()+
  geom_smooth(method=lm , se=TRUE, formula = y ~ x,col="black") +
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")

pl_lin_vss <- ggplot(allo_sm_trans_nca_iv_sd,aes(x=BW,y=vss.pred))+
  geom_point(aes(col=as.factor(DOSE_MGKG))) +
  scale_x_log10()+scale_y_log10()+
  geom_smooth(method=lm , se=TRUE, formula = y ~ x,col="black") +
  labs(x="Body weight (kg)",col="Dose (mg/kg)",subtitle=paste0("Slope = Allometric coefficient = ",ci.vss))+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_lin_vss

pl_lin_halflife <- ggplot(allo_sm_trans_nca_iv_sd,aes(x=BW,y=half.life))+
  geom_point(aes(col=as.factor(DOSE_MGKG))) +
  scale_x_log10()+scale_y_log10()+
  geom_smooth(method=lm , se=TRUE, formula = y ~ x,col="black") +
  labs(x="Body weight (kg)",col="Dose (mg/kg)",subtitle=paste0("Slope = Allometric coefficient = ",ci.half.life))+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_lin_halflife


# Step 2: results - iv po -------------------------------------------------

pl_rat_pk_ev <- ggplot(subset(mutate(trans_sd_err,ROUTE=factor(ifelse(ROUTE==1,"IV","PO"))),C1_BQL==0&CMPD=="SM"&SPEC%in%c("Rat")),
                             aes(x=time,y=C1_Y,col=as.factor(ROUTE)))+
  geom_point() +
  # geom_line(aes(group=ID)) +
  stat_summary(fun=median, geom='line', ) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,168,24))+
  facet_wrap(~DOSE_MGKG,labeller="label_both")+
  geom_hline(yintercept = 0.05)+
  labs(col="Route",subtitle="iv/po sd SM rat",x="Time (h)",y="Cplasma (?g/mL)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[c(2,1)])+ theme(legend.position="bottom")
pl_rat_pk_ev

pl_rat_pk_ev_zoom <- ggplot(subset(mutate(trans_sd_err,ROUTE=factor(ifelse(ROUTE==1,"IV","PO"))),C1_BQL==0&CMPD=="SM"&SPEC%in%c("Rat")),
                       aes(x=time,y=C1_Y,col=as.factor(ROUTE)))+
  geom_point() +
  geom_line(aes(group=ID)) +
  # stat_summary(fun=median, geom='line', ) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,168,8))+
  coord_cartesian(xlim=c(0,48))+
  facet_wrap(~DOSE_MGKG,labeller="label_both")+
  geom_hline(yintercept = 0.05)+
  labs(col="Route",subtitle="iv/po sd SM rat",x="Time (h)",y="Cplasma (?g/mL)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[c(2,1)])+ theme(legend.position="bottom")
pl_rat_pk_ev_zoom

pl_rat_nca_ev <- ggplot(subset(mutate(nca_sd_trans,ROUTE=factor(ifelse(ROUTE==1,"IV","PO"))),CMPD=="SM"&SPEC=="Rat"),
       aes(x=DOSE_MGKG,y=cl.pred*1000,col=ROUTE))+
  geom_point() +
  geom_smooth(method=lm , se=TRUE, formula = y ~ x) +
  labs(col="Route",subtitle="iv/po sd SM rat",x="Dose (mg/kg)",y="CL (mL/h)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_rat_nca_ev

# Step 3: results - mAb nonlinearity --------------------------------------

pl_trans_pk_nonlin <- ggplot(subset(trans_sd_err,C1_BQL==0&CMPD=="mAb"&ROUTE==1&SPEC%in%c("Rat","Monkey")),
       aes(x=time/24,y=C1_Y,col=as.factor(ADA)))+
  geom_point() +
  geom_line(aes(group=ID)) +
  stat_summary(fun=median, geom='line', ) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  facet_wrap(~SPEC+DOSE_MGKG,labeller="label_both")+
  geom_hline(yintercept = 0.05)+
  labs(col="ADA",subtitle="iv sd mAb",x="Time (days)",y="Cplasma (?g/mL)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[c(2,1)])+ theme(legend.position="bottom")
pl_trans_pk_nonlin

pl_trans_nca_nonlin <- ggplot(subset(nca_sd_trans,CMPD=="mAb"&ROUTE==1&SPEC%in%c("Rat","Monkey")),
                              aes(x=DOSE_MGKG,y=cl.pred*1000,col=as.factor(ADA)))+
  geom_point() +
  stat_summary(fun=median, geom='line', ) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  facet_wrap(~SPEC)+
  labs(col="ADA",subtitle="iv sd mAb",x="Dose (mg/kg)",y="CL (mL/h)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[c(2,1)])+ theme(legend.position="bottom")
pl_trans_nca_nonlin

# Step 4: PKPD - SELECT TPs ------------------------------
trans_pkpd_q2d_err <- read.csv("outsel_err_n10.csv") %>%
  subset(SPEC%in%c("Monkey","Rat")&NDOS>1&ROUTE==2)

# Time point selection per variable
trans_pkpd_q2d_c1 <- dplyr::select(trans_pkpd_q2d_err,ID,TIME=time,DV=C1_Y,BQL=C1_BQL) %>%
  subset(TIME%in%c(0.5,4,24,24*seq(2,14,4))) %>% mutate(TYPE="Cplasma",BSL=0)

trans_pkpd_q2d_cb <- dplyr::select(trans_pkpd_q2d_err,ID,TIME=time,DV=CB_Y,BQL=CB_BQL) %>%
  subset(TIME%in%c(14*24)) %>% mutate(TYPE="Cbrain",BSL=0) # terminal

trans_pkpd_q2d_eff <- dplyr::select(trans_pkpd_q2d_err,ID,TIME=time,DV=EFF_Y,BQL=EFF_BQL) %>%
  subset(TIME%in%c(4,316,seq(0,14*24,2*24))) %>% mutate(TYPE="Efficacy_Biomarker")
trans_pkpd_q2d_eff <- dplyr::left_join(trans_pkpd_q2d_eff,dplyr::select(subset(trans_pkpd_q2d_eff,TIME==0),ID,BSL=DV))

trans_pkpd_q2d_saf <- dplyr::select(trans_pkpd_q2d_err,ID,TIME=time,DV=SAF_Y,BQL=SAF_BQL) %>%
  subset(TIME%in%c(4,316,seq(0,14*24,2*24))) %>% mutate(TYPE="Safety_Biomarker")
trans_pkpd_q2d_saf <- dplyr::left_join(trans_pkpd_q2d_saf,dplyr::select(subset(trans_pkpd_q2d_saf,TIME==0),ID,BSL=DV))

trans_pkpd_q2d_long <- rbind(trans_pkpd_q2d_c1,trans_pkpd_q2d_cb,trans_pkpd_q2d_eff,trans_pkpd_q2d_saf) %>%
  dplyr::left_join(dplyr::select(subset(trans_pkpd_q2d_err,!duplicated(ID)),ID,SPEC,DOSE_MGKG,ROUTE,CMPD,BW,ADA))

trans_pkpd_q2d_wide <- subset(trans_pkpd_q2d_long,BQL==0) %>% 
  mutate(BQL=NULL,BSL=NULL) %>%
  tidyr::spread(key=TYPE, value=DV)

trans_pkpd_q2d_wide0 <- subset(trans_pkpd_q2d_wide,TIME==0) %>%
  mutate(Efficacy_Baseline=Efficacy_Biomarker,Safety_Baseline=Safety_Biomarker) %>%
  dplyr::select(ID,Efficacy_Baseline,Safety_Baseline)

trans_pkpd_q2d_wide <- dplyr::left_join(trans_pkpd_q2d_wide,trans_pkpd_q2d_wide0) %>%
  mutate(Efficacy_Percchange=100*(Efficacy_Biomarker-Efficacy_Baseline)/Efficacy_Baseline,
         Safety_Percchange=100*(Safety_Biomarker-Safety_Baseline)/Safety_Baseline)

write.csv(trans_pkpd_q2d_long, file="./HO/Step4/trans_pkpd_q2d_long.csv", row.names = FALSE, quote = FALSE)
write.csv(trans_pkpd_q2d_wide, file="./HO/Step4/trans_pkpd_q2d_wide.csv", row.names = FALSE, quote = FALSE)

# Step 4: Results - PKPD --------------------------------------------------

pl_trans_pkpd_pk <- ggplot(subset(trans_pkpd_q2d_long,DV>0.05&TYPE%in%c("Cplasma","Cbrain")))+
  geom_point(aes(x=TIME/24,y=DV,group=ID,col=CMPD),alpha=0.5)+
  geom_line(aes(x=TIME/24,y=DV,group=ID,col=CMPD),alpha=0.5)+
  scale_y_log10()+
  facet_wrap(~SPEC+DOSE_MGKG+TYPE,labeller = "label_both")+
  geom_hline(yintercept = 0.05)+
  scale_x_continuous(breaks=seq(0,100,2))+
  labs(x="Time (days)",y="Concentration (?g/mL)")+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_trans_pkpd_pk

pl_trans_pkpd_eff <- ggplot(subset(trans_pkpd_q2d_long,TYPE%in%c("Efficacy_Biomarker")))+
  geom_point(aes(x=TIME/24,y=DV,group=ID,col=CMPD),alpha=0.5)+
  geom_line(aes(x=TIME/24,y=DV,group=ID,col=CMPD),alpha=0.5)+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_x_continuous(breaks=seq(0,100,2))+
  labs(x="Time (days)",y="Biomarker (absolute value)",subtitle = 'Efficacy biomarker')+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_trans_pkpd_eff

pl_trans_pkpd_eff_chg <- ggplot(subset(trans_pkpd_q2d_long,TYPE%in%c("Efficacy_Biomarker")))+
  geom_point(aes(x=TIME/24,y=100*(DV-BSL)/BSL,group=ID,col=CMPD),alpha=0.5)+
  geom_line(aes(x=TIME/24,y=100*(DV-BSL)/BSL,group=ID,col=CMPD),alpha=0.5)+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_x_continuous(breaks=seq(0,100,2))+
  labs(x="Time (days)",y="Biomarker (% change from baseline)",subtitle = 'Efficacy biomarker')+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_trans_pkpd_eff_chg

pl_trans_pkpd_saf <- ggplot(subset(trans_pkpd_q2d_long,TYPE%in%c("Safety_Biomarker")))+
  geom_point(aes(x=TIME/24,y=DV,group=ID,col=CMPD),alpha=0.5)+
  geom_line(aes(x=TIME/24,y=DV,group=ID,col=CMPD),alpha=0.5)+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_x_continuous(breaks=seq(0,100,2))+
  labs(x="Time (days)",y="Biomarker (absolute value)",subtitle = 'Safety biomarker')+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_trans_pkpd_saf

pl_trans_pkpd_saf_chg <- ggplot(subset(trans_pkpd_q2d_long,TYPE%in%c("Safety_Biomarker")))+
  geom_point(aes(x=TIME/24,y=100*(DV-BSL)/BSL,group=ID,col=CMPD),alpha=0.5)+
  geom_line(aes(x=TIME/24,y=100*(DV-BSL)/BSL,group=ID,col=CMPD),alpha=0.5)+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  scale_x_continuous(breaks=seq(0,100,2))+
  labs(x="Time (days)",y="Biomarker (% change from baseline)",subtitle = 'Safety biomarker')+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_trans_pkpd_saf_chg

pl_trans_pkpd_window_t336 <- ggplot(subset(trans_pkpd_q2d_wide,TIME==336))+
  geom_point(aes(x=Efficacy_Percchange,y=Safety_Percchange,group=ID,col=CMPD),alpha=0.5)+
  # geom_line(aes(x=Efficacy_Percchange,y=Safety_Percchange,group=ID,col=CMPD))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  coord_cartesian(xlim=c(-100,40),ylim=c(-100,40))+
  scale_x_continuous(breaks=seq(-100,100,20))+
  scale_y_continuous(breaks=seq(-100,100,20))+
  geom_hline(yintercept = -20,col="red")+
  geom_vline(xintercept = -40,col="darkgreen")+
  labs(subtitle=paste0("subset at time = ",unique(subset(trans_pkpd_q2d_wide,TIME==336)$TIME)/24," days"))+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_trans_pkpd_window_t336


pl_trans_pkpd_window_t316 <- ggplot(subset(trans_pkpd_q2d_wide,TIME==316))+
  geom_point(aes(x=Efficacy_Percchange,y=Safety_Percchange,group=ID,col=CMPD),alpha=0.5)+
  # geom_line(aes(x=Efficacy_Percchange,y=Safety_Percchange,group=ID,col=CMPD))+
  facet_wrap(~SPEC+DOSE_MGKG,labeller = "label_both")+
  coord_cartesian(xlim=c(-100,40),ylim=c(-100,40))+
  scale_x_continuous(breaks=seq(-100,100,20))+
  scale_y_continuous(breaks=seq(-100,100,20))+
  geom_hline(yintercept = -20,col="red")+
  geom_vline(xintercept = -40,col="darkgreen")+
  labs(subtitle=paste0("subset at time = ",signif(unique(subset(trans_pkpd_q2d_wide,TIME==316)$TIME)/24,4)," days"))+
  theme_lapp()+scale_color_manual(values=col_lapp()[])+ theme(legend.position="bottom")
pl_trans_pkpd_window_t316


# Step 5: human POC simulations PKPD ------------------------------

# User input

nid <- 200

s_df <- expand.grid(SPEC = c(1),
                    DOSE_MGKG = c(5,20,50),
                    ROUTE = c(2),
                    # II = c(48),
                    # NDOS = c(1,10),
                    CMPD = c(1,2)
) %>%
  mutate(
    II=ifelse(CMPD==1,7*24,24),
    DOSE_MGKG=ifelse(CMPD==1,DOSE_MGKG,DOSE_MGKG/10),
    NDOS=ifelse(CMPD==1,4,28),
    TVCL=ifelse(CMPD==1,p_mab[["TVCL"]],p_sm[["TVCL"]]),
    TVV1=ifelse(CMPD==1,p_mab[["TVV1"]],p_sm[["TVV1"]]),
    TVKA=ifelse(CMPD==1,p_mab[["TVKA"]],p_sm[["TVKA"]]),
    TVF=ifelse(CMPD==1,p_mab[["TVF"]],p_sm[["TVF"]]),
    TVQ=ifelse(CMPD==1,p_mab[["TVQ"]],p_sm[["TVQ"]]),
    TVV2=ifelse(CMPD==1,p_mab[["TVV2"]],p_sm[["TVV2"]]),
    TVVMAX=ifelse(CMPD==1,p_mab[["TVVMAX"]],p_sm[["TVVMAX"]]),
    TVVMAX=ifelse(CMPD==1&SPEC%in%c(3,4),0*TVVMAX,TVVMAX),        # Rodent mab non linearity not picked up
    TVKM=ifelse(CMPD==1,p_mab[["TVKM"]],p_sm[["TVKM"]]),
    KP=ifelse(CMPD==1,p_mab[["KP"]],p_sm[["KP"]]),
    TVQB=ifelse(CMPD==1,p_mab[["TVQB"]],p_sm[["TVQB"]]),
    TVVB=ifelse(CMPD==1,p_mab[["TVVB"]],p_sm[["TVVB"]]),
    TV_SAF_KOUT=TV_SAF_KOUT,TV_SAF_KIN=TV_SAF_KIN,TV_EFF_KOUT=TV_EFF_KOUT,TV_EFF_KIN=TV_EFF_KIN,
    SAF_MAX=ifelse(CMPD==1,p_mab[["SAF_MAX"]],p_sm[["SAF_MAX"]]),
    SAF_C50=ifelse(CMPD==1,p_mab[["SAF_C50"]],p_sm[["SAF_C50"]]),
    SAF_GAM=ifelse(CMPD==1,p_mab[["SAF_GAM"]],p_sm[["SAF_GAM"]]),
    EFF_MAX=ifelse(CMPD==1,p_mab[["EFF_MAX"]],p_sm[["EFF_MAX"]]),
    EFF_C50=ifelse(CMPD==1,p_mab[["EFF_C50"]],p_sm[["EFF_C50"]]),
    EFF_C50=ifelse(CMPD==1&SPEC%in%c(3,4),10*EFF_C50,EFF_C50),   # Rodent EC50 higher for mab
    EFF_GAM=ifelse(CMPD==1,p_mab[["EFF_GAM"]],p_sm[["EFF_GAM"]]),
    BW=ifelse(SPEC==1,70,ifelse(SPEC==2,3,ifelse(SPEC==3,0.3,ifelse(SPEC==4,0.02,-999)))), # -999 to ensure error as not wanted input
    DOSE_MG=DOSE_MGKG*BW)
s_df$SCEN <- 1:nrow(s_df)
nrep        <- nrow(s_df)

# eta generation
samp_shell <- data.frame(ID= 1:(nid*nrep), SCEN = rep(1:nrep,nid), REP = rep(1:nid,each=nrep))
set.seed(1) ; eta1 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var1))
set.seed(2) ; eta2 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var2))
set.seed(3) ; eta3 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var3))
set.seed(4) ; eta4 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var4))
set.seed(5) ; eta5 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var5))
set.seed(6) ; eta6 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var6))
samp_ETA <-  mutate(samp_shell,ETA1=eta1,ETA2=eta2,ETA3=eta3,ETA4=eta4,ETA5=eta5,ETA6=eta6)

set.seed(10) ; adadraw <- runif(nrow(samp_ETA))
samp <- dplyr::left_join(samp_ETA,s_df,by="SCEN") %>%
  mutate(ADADRAW=adadraw,ADA=ifelse(CMPD==1&SPEC==2&ADADRAW<prob_monkey_ada,1,ifelse(CMPD==1&SPEC%in%c(3,4)&ADADRAW<prob_rodent_ada,1,0)))

# Testing via # x <- 1 # unlist(samp[x,])
out <- lapply(1:nrow(samp), function(x){
  print(x)
  parm  <- par_func(TVCL=samp[x,]$TVCL,TVV1=samp[x,]$TVV1,TVQ=samp[x,]$TVQ,TVV2=samp[x,]$TVV2,TVVMAX=samp[x,]$TVVMAX,TVKM=samp[x,]$TVKM,TVKA=samp[x,]$TVKA,
                    ETA1=samp[x,]$ETA1,ETA2=samp[x,]$ETA2,ETA3=samp[x,]$ETA3,ETA4=samp[x,]$ETA4,ETA5=samp[x,]$ETA5,ETA6=samp[x,]$ETA6,
                    KP=samp[x,]$KP,TVQB=samp[x,]$TVQB,TVVB=samp[x,]$TVVB,
                    TV_SAF_KOUT=samp[x,]$TV_SAF_KOUT,TV_SAF_KIN=samp[x,]$TV_SAF_KIN,TV_EFF_KOUT=samp[x,]$TV_EFF_KOUT,TV_EFF_KIN=samp[x,]$TV_EFF_KIN,
                    SAF_MAX=samp[x,]$SAF_MAX,SAF_C50=samp[x,]$SAF_C50,SAF_GAM=samp[x,]$SAF_GAM,
                    EFF_MAX=samp[x,]$EFF_MAX,EFF_C50=samp[x,]$EFF_C50,EFF_GAM=samp[x,]$EFF_GAM,
                    SPEC=samp[x,]$SPEC,
                    BW=samp[x,]$BW,
                    ADA=samp[x,]$ADA,
                    ADASTART=ada_onset) # If ADA=1, ADA effect kicks in onset of time, currently no IIV
  inits <- init_func(parm)
  if(unlist(samp[x,])[["ROUTE"]]==1) evnt  <- dose_func(cmt=1,
                                                        value=samp[x,]$DOSE_MG,
                                                        # tinf=2, # careful for infusion dose in cmt=3, for bolus dose in cmt=1
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  if(unlist(samp[x,])[["ROUTE"]]!=1) evnt  <- dose_func(cmt=5,
                                                        value=samp[x,]$DOSE_MG * samp[x,]$TVF,
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  
  times  <- sort(c(1,seq(0,4*7*24,6)))
  data.frame(deSolve::lsoda(inits,times,des_func,parm,events=list(data=evnt)),
             ID=unlist(samp[x,])[["ID"]]
  )
})
out <- do.call(rbind,out)

outsel_human <- dplyr::select(out,time,C1,AUC=A4,CB,SAFETY=A7,EFFICACY=A8,ID) %>% 
  dplyr::left_join(dplyr::select(samp,ID:REP,SPEC:CMPD,BW,DOSE_MG,ADA),by="ID") %>%
  mutate(SPEC=ifelse(SPEC==1,"Human",ifelse(SPEC==2,"Monkey",ifelse(SPEC==3,"Rat",ifelse(SPEC==4,"Mouse","ERROR")))),
         CMPD=as.factor(ifelse(CMPD==1,"mAb","SM")))

outsel_human0 <- subset(outsel_human,time==0) %>%
  mutate(SAFETY0 = SAFETY, EFFICACY0 = EFFICACY) %>%
  dplyr::select(ID,SAFETY0,EFFICACY0)

outsel_human <- dplyr::left_join(outsel_human,outsel_human0) %>% 
  mutate(EFF_CHG = 100 * (EFFICACY-EFFICACY0)/EFFICACY0,
         SAF_CHG = 100 * (SAFETY-SAFETY0)/SAFETY0)

write.csv(outsel_human, file=paste0("outsel_human","_n",nid,".csv"), row.names = FALSE, quote = FALSE, na = ".")

sum_human <- outsel_human %>% 
  group_by(CMPD,DOSE_MGKG,time) %>% 
  dplyr::summarise(N=length(C1), 
                   C1_05 = quantile(C1,probs=c(0.05)),
                   C1_50 = quantile(C1,probs=c(0.5)),
                   C1_95 = quantile(C1,probs=c(0.95)),
                   CB_05 = quantile(CB,probs=c(0.05)),
                   CB_50 = quantile(CB,probs=c(0.5)),
                   CB_95 = quantile(CB,probs=c(0.95)),
                   EFF_CHG_05 = quantile(EFF_CHG,probs=c(0.05)),
                   EFF_CHG_50 = quantile(EFF_CHG,probs=c(0.5)),
                   EFF_CHG_95 = quantile(EFF_CHG,probs=c(0.95)),
                   SAF_CHG_05 = quantile(SAF_CHG,probs=c(0.05)),
                   SAF_CHG_50 = quantile(SAF_CHG,probs=c(0.5)),
                   SAF_CHG_95 = quantile(SAF_CHG,probs=c(0.95))
                   ) %>% as.data.frame()

pl_hum_cp <- ggplot(data=subset(sum_human))+
  facet_wrap(~CMPD)+
  geom_ribbon(aes(x= time/24, ymin = C1_05, ymax = C1_95,fill=as.factor(DOSE_MGKG)),alpha=0.5) +
  geom_line(aes(x=time/24,y=C1_50,col=as.factor(DOSE_MGKG)))+
  scale_y_log10() +  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,100,4))+
  geom_hline(yintercept=0.05,col="grey50")+
  coord_cartesian(ylim=c(0.05,1000))+
  scale_color_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  scale_fill_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  labs(col="Dose (mg/kg)",fill="Dose (mg/kg)", x="Time (days)",y="Cplasma (?g/mL)",subtitle="Human predictions: mAb sc qw dosing; SM po qd dosing")+
  theme_lapp()+theme(legend.position="bottom")
pl_hum_cp

pl_hum_cb <- ggplot(data=subset(sum_human))+
  facet_wrap(~CMPD)+
  geom_ribbon(aes(x= time/24, ymin = CB_05, ymax = CB_95,fill=as.factor(DOSE_MGKG)),alpha=0.5) +
  geom_line(aes(x=time/24,y=CB_50,col=as.factor(DOSE_MGKG)))+
  scale_y_log10() +  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,100,4))+
  geom_hline(yintercept=0.05,col="grey50")+
  coord_cartesian(ylim=c(0.05,1000))+
  scale_color_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  scale_fill_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  labs(col="Dose (mg/kg)",fill="Dose (mg/kg)", x="Time (days)",y="Cbrain (?g/mL)",subtitle="Human predictions: mAb sc qw dosing; SM po qd dosing")+
  theme_lapp()+theme(legend.position="bottom")
pl_hum_cb

pl_hum_eff <- ggplot(data=subset(sum_human))+
  facet_wrap(~CMPD)+
  geom_ribbon(aes(x= time/24, ymin = EFF_CHG_05, ymax = EFF_CHG_95,fill=as.factor(DOSE_MGKG)),alpha=0.5) +
  geom_line(aes(x=time/24,y=EFF_CHG_50,col=as.factor(DOSE_MGKG)))+
  scale_y_continuous(breaks=seq(-100,0,10))+
  scale_x_continuous(breaks=seq(0,100,4))+
  coord_cartesian(ylim=c(-100,0))+
  geom_hline(yintercept=0.05,col="grey50")+
  scale_color_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  scale_fill_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  labs(col="Dose (mg/kg)",fill="Dose (mg/kg)", x="Time (days)",y="Efficacy biomarker (% change from baseline)",subtitle="Human predictions: mAb sc qw dosing; SM po qd dosing")+
  theme_lapp()+theme(legend.position="bottom")
  pl_hum_eff

pl_hum_saf <- ggplot(data=subset(sum_human))+
  facet_wrap(~CMPD)+
  geom_ribbon(aes(x= time/24, ymin = SAF_CHG_05, ymax = SAF_CHG_95,fill=as.factor(DOSE_MGKG)),alpha=0.5) +
  geom_line(aes(x=time/24,y=SAF_CHG_50,col=as.factor(DOSE_MGKG)))+
  scale_y_continuous(breaks=seq(-100,0,10))+
  scale_x_continuous(breaks=seq(0,100,4))+
  coord_cartesian(ylim=c(-100,0))+
  geom_hline(yintercept=0.05,col="grey50")+
  scale_color_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  scale_fill_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  labs(col="Dose (mg/kg)",fill="Dose (mg/kg)", x="Time (days)",y="Safety biomarker (% change from baseline)",subtitle="Human predictions: mAb sc qw dosing; SM po qd dosing")+
  theme_lapp()+theme(legend.position="bottom")
pl_hum_saf


# Step 6: Typical Cyno AUC at NOAEL ------------------------------

# User input - Note: ada set to 0, eta's set to 0, nid set to 1, 28w dosing + period (NDOS, and times)

nid <- 1

s_df <- expand.grid(SPEC = c(2),
                    DOSE_MGKG = c(100),
                    ROUTE = c(2),
                    NDOS = c(1,28),
                    CMPD = c(1)
) %>%
  mutate(
    II=ifelse(CMPD==1,7*24,24),
    DOSE_MGKG=ifelse(CMPD==1,DOSE_MGKG,DOSE_MGKG/10),
    # NDOS=ifelse(CMPD==1,4,28),
    TVCL=ifelse(CMPD==1,p_mab[["TVCL"]],p_sm[["TVCL"]]),
    TVV1=ifelse(CMPD==1,p_mab[["TVV1"]],p_sm[["TVV1"]]),
    TVKA=ifelse(CMPD==1,p_mab[["TVKA"]],p_sm[["TVKA"]]),
    TVF=ifelse(CMPD==1,p_mab[["TVF"]],p_sm[["TVF"]]),
    TVQ=ifelse(CMPD==1,p_mab[["TVQ"]],p_sm[["TVQ"]]),
    TVV2=ifelse(CMPD==1,p_mab[["TVV2"]],p_sm[["TVV2"]]),
    TVVMAX=ifelse(CMPD==1,p_mab[["TVVMAX"]],p_sm[["TVVMAX"]]),
    TVVMAX=ifelse(CMPD==1&SPEC%in%c(3,4),0*TVVMAX,TVVMAX),        # Rodent mab non linearity not picked up
    TVKM=ifelse(CMPD==1,p_mab[["TVKM"]],p_sm[["TVKM"]]),
    KP=ifelse(CMPD==1,p_mab[["KP"]],p_sm[["KP"]]),
    TVQB=ifelse(CMPD==1,p_mab[["TVQB"]],p_sm[["TVQB"]]),
    TVVB=ifelse(CMPD==1,p_mab[["TVVB"]],p_sm[["TVVB"]]),
    TV_SAF_KOUT=TV_SAF_KOUT,TV_SAF_KIN=TV_SAF_KIN,TV_EFF_KOUT=TV_EFF_KOUT,TV_EFF_KIN=TV_EFF_KIN,
    SAF_MAX=ifelse(CMPD==1,p_mab[["SAF_MAX"]],p_sm[["SAF_MAX"]]),
    SAF_C50=ifelse(CMPD==1,p_mab[["SAF_C50"]],p_sm[["SAF_C50"]]),
    SAF_GAM=ifelse(CMPD==1,p_mab[["SAF_GAM"]],p_sm[["SAF_GAM"]]),
    EFF_MAX=ifelse(CMPD==1,p_mab[["EFF_MAX"]],p_sm[["EFF_MAX"]]),
    EFF_C50=ifelse(CMPD==1,p_mab[["EFF_C50"]],p_sm[["EFF_C50"]]),
    EFF_C50=ifelse(CMPD==1&SPEC%in%c(3,4),10*EFF_C50,EFF_C50),   # Rodent EC50 higher for mab
    EFF_GAM=ifelse(CMPD==1,p_mab[["EFF_GAM"]],p_sm[["EFF_GAM"]]),
    BW=ifelse(SPEC==1,70,ifelse(SPEC==2,3,ifelse(SPEC==3,0.3,ifelse(SPEC==4,0.02,-999)))), # -999 to ensure error as not wanted input
    DOSE_MG=DOSE_MGKG*BW)
s_df$SCEN <- 1:nrow(s_df)
nrep        <- nrow(s_df)

# eta generation: set to 0 for Typical sim
samp_shell <- data.frame(ID= 1:(nid*nrep), SCEN = rep(1:nrep,nid), REP = rep(1:nid,each=nrep))
set.seed(1) ; eta1 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var1))
set.seed(2) ; eta2 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var2))
set.seed(3) ; eta3 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var3))
set.seed(4) ; eta4 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var4))
set.seed(5) ; eta5 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var5))
set.seed(6) ; eta6 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var6))
# samp_ETA <-  mutate(samp_shell,ETA1=eta1,ETA2=eta2,ETA3=eta3,ETA4=eta4,ETA5=eta5,ETA6=eta6)
samp_ETA <-  mutate(samp_shell,ETA1=0,ETA2=0,ETA3=0,ETA4=0,ETA5=0,ETA6=0)

set.seed(10) ; adadraw <- runif(nrow(samp_ETA))
samp <- dplyr::left_join(samp_ETA,s_df,by="SCEN") %>%
  # mutate(ADADRAW=adadraw,ADA=ifelse(CMPD==1&SPEC==2&ADADRAW<prob_monkey_ada,1,ifelse(CMPD==1&SPEC%in%c(3,4)&ADADRAW<prob_rodent_ada,1,0)))
  mutate(ADADRAW=1,ADA=ifelse(CMPD==1&SPEC==2&ADADRAW<prob_monkey_ada,1,ifelse(CMPD==1&SPEC%in%c(3,4)&ADADRAW<prob_rodent_ada,1,0)))

# Testing via # x <- 1 # unlist(samp[x,])
out <- lapply(1:nrow(samp), function(x){
  print(x)
  parm  <- par_func(TVCL=samp[x,]$TVCL,TVV1=samp[x,]$TVV1,TVQ=samp[x,]$TVQ,TVV2=samp[x,]$TVV2,TVVMAX=samp[x,]$TVVMAX,TVKM=samp[x,]$TVKM,TVKA=samp[x,]$TVKA,
                    ETA1=samp[x,]$ETA1,ETA2=samp[x,]$ETA2,ETA3=samp[x,]$ETA3,ETA4=samp[x,]$ETA4,ETA5=samp[x,]$ETA5,ETA6=samp[x,]$ETA6,
                    KP=samp[x,]$KP,TVQB=samp[x,]$TVQB,TVVB=samp[x,]$TVVB,
                    TV_SAF_KOUT=samp[x,]$TV_SAF_KOUT,TV_SAF_KIN=samp[x,]$TV_SAF_KIN,TV_EFF_KOUT=samp[x,]$TV_EFF_KOUT,TV_EFF_KIN=samp[x,]$TV_EFF_KIN,
                    SAF_MAX=samp[x,]$SAF_MAX,SAF_C50=samp[x,]$SAF_C50,SAF_GAM=samp[x,]$SAF_GAM,
                    EFF_MAX=samp[x,]$EFF_MAX,EFF_C50=samp[x,]$EFF_C50,EFF_GAM=samp[x,]$EFF_GAM,
                    SPEC=samp[x,]$SPEC,
                    BW=samp[x,]$BW,
                    ADA=samp[x,]$ADA,
                    ADASTART=ada_onset) # If ADA=1, ADA effect kicks in onset of time, currently no IIV
  inits <- init_func(parm)
  if(unlist(samp[x,])[["ROUTE"]]==1) evnt  <- dose_func(cmt=1,
                                                        value=samp[x,]$DOSE_MG,
                                                        # tinf=2, # careful for infusion dose in cmt=3, for bolus dose in cmt=1
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  if(unlist(samp[x,])[["ROUTE"]]!=1) evnt  <- dose_func(cmt=5,
                                                        value=samp[x,]$DOSE_MG * samp[x,]$TVF,
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  
  times  <- sort(c(1,seq(0,29*7*24,1)))
  # times  <- c(0,1*7*24,27*7*24,28*7*24)
  
  data.frame(deSolve::lsoda(inits,times,des_func,parm,events=list(data=evnt)),
             ID=unlist(samp[x,])[["ID"]]
  )
})
out <- do.call(rbind,out)

outsel_noael <- dplyr::select(out,time,C1,AUC=A4,CB,SAFETY=A7,EFFICACY=A8,ID) %>% 
  dplyr::left_join(dplyr::select(samp,ID:REP,SPEC:CMPD,BW,DOSE_MG,ADA),by="ID") %>%
  mutate(SPEC=ifelse(SPEC==1,"Human",ifelse(SPEC==2,"Monkey",ifelse(SPEC==3,"Rat",ifelse(SPEC==4,"Mouse","ERROR")))),
         CMPD=as.factor(ifelse(CMPD==1,"mAb","SM")))

ggplot(data=subset(outsel_noael,NDOS==1))+
  facet_wrap(~CMPD)+
  geom_line(aes(x=time/24/7,y=C1,col=as.factor(DOSE_MGKG)))+
  scale_y_log10() +  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,100,2))+
  coord_cartesian(ylim=c(100,10000))+
  geom_hline(yintercept=0.05,col="grey50")+
  scale_color_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  scale_fill_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  labs(col="Dose (mg/kg)",fill="Dose (mg/kg)", x="Time (days)",y="Cplasma (?g/mL)")+
  theme_lapp()+theme(legend.position="bottom")

ggplot(data=subset(outsel_noael,NDOS!=1))+
  facet_wrap(~CMPD)+
  geom_line(aes(x=time/24/7,y=C1,col=as.factor(DOSE_MGKG)))+
  scale_y_log10() +  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks=seq(0,100,2))+
  coord_cartesian(ylim=c(100,10000))+
  geom_hline(yintercept=0.05,col="grey50")+
  scale_color_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  scale_fill_manual(values=col_lapp()[c(1,2,3,6,7,8)])+
  labs(col="Dose (mg/kg)",fill="Dose (mg/kg)", x="Time (days)",y="Cplasma (?g/mL)")+
  theme_lapp()+theme(legend.position="bottom")

auc_noael_sd <- subset(outsel_noael,NDOS==1&time==28*7*24)$AUC
signif(subset(outsel_noael,NDOS==1&time==28*7*24)$AUC - subset(outsel_noael,NDOS==1&time==27*7*24)$AUC,2)
# check: AUC0-27d =AUC0-28d so can be named AUC0-inf

auc_noael_w1 <- subset(outsel_noael,NDOS!=1&time==1*7*24)$AUC
auc_noael_ss <- subset(outsel_noael,NDOS!=1&time==28*7*24)$AUC - subset(outsel_noael,NDOS!=1&time==27*7*24)$AUC

auc_noael_sd ; auc_noael_w1 ; auc_noael_ss
# target AUC with margin
auc_noael_sd/10 ; auc_noael_w1/10 ; auc_noael_ss/10




# Step 7 : AUC based HED and MABEL ------------------------------

# User input - Note: ada set to 0, eta's set to 0, nid set to 1, 28w dosing + period (NDOS, and times)

nid <- 1

s_df <- expand.grid(SPEC = c(1),
                    DOSE_MGKG = seq(1,150,1)/10,
                    ROUTE = c(2),
                    NDOS = c(1,28),
                    CMPD = c(1)
) %>%
  mutate(
    II=ifelse(CMPD==1,7*24,24),
    DOSE_MGKG=ifelse(CMPD==1,DOSE_MGKG,DOSE_MGKG/10),
    # NDOS=ifelse(CMPD==1,4,28),
    TVCL=ifelse(CMPD==1,p_mab[["TVCL"]],p_sm[["TVCL"]]),
    TVV1=ifelse(CMPD==1,p_mab[["TVV1"]],p_sm[["TVV1"]]),
    TVKA=ifelse(CMPD==1,p_mab[["TVKA"]],p_sm[["TVKA"]]),
    TVF=ifelse(CMPD==1,p_mab[["TVF"]],p_sm[["TVF"]]),
    TVQ=ifelse(CMPD==1,p_mab[["TVQ"]],p_sm[["TVQ"]]),
    TVV2=ifelse(CMPD==1,p_mab[["TVV2"]],p_sm[["TVV2"]]),
    TVVMAX=ifelse(CMPD==1,p_mab[["TVVMAX"]],p_sm[["TVVMAX"]]),
    TVVMAX=ifelse(CMPD==1&SPEC%in%c(3,4),0*TVVMAX,TVVMAX),        # Rodent mab non linearity not picked up
    TVKM=ifelse(CMPD==1,p_mab[["TVKM"]],p_sm[["TVKM"]]),
    KP=ifelse(CMPD==1,p_mab[["KP"]],p_sm[["KP"]]),
    TVQB=ifelse(CMPD==1,p_mab[["TVQB"]],p_sm[["TVQB"]]),
    TVVB=ifelse(CMPD==1,p_mab[["TVVB"]],p_sm[["TVVB"]]),
    TV_SAF_KOUT=TV_SAF_KOUT,TV_SAF_KIN=TV_SAF_KIN,TV_EFF_KOUT=TV_EFF_KOUT,TV_EFF_KIN=TV_EFF_KIN,
    SAF_MAX=ifelse(CMPD==1,p_mab[["SAF_MAX"]],p_sm[["SAF_MAX"]]),
    SAF_C50=ifelse(CMPD==1,p_mab[["SAF_C50"]],p_sm[["SAF_C50"]]),
    SAF_GAM=ifelse(CMPD==1,p_mab[["SAF_GAM"]],p_sm[["SAF_GAM"]]),
    EFF_MAX=ifelse(CMPD==1,p_mab[["EFF_MAX"]],p_sm[["EFF_MAX"]]),
    EFF_C50=ifelse(CMPD==1,p_mab[["EFF_C50"]],p_sm[["EFF_C50"]]),
    EFF_C50=ifelse(CMPD==1&SPEC%in%c(3,4),10*EFF_C50,EFF_C50),   # Rodent EC50 higher for mab
    EFF_GAM=ifelse(CMPD==1,p_mab[["EFF_GAM"]],p_sm[["EFF_GAM"]]),
    BW=ifelse(SPEC==1,70,ifelse(SPEC==2,3,ifelse(SPEC==3,0.3,ifelse(SPEC==4,0.02,-999)))), # -999 to ensure error as not wanted input
    DOSE_MG=DOSE_MGKG*BW)
s_df$SCEN <- 1:nrow(s_df)
nrep        <- nrow(s_df)

# eta generation: set to 0 for Typical sim
samp_shell <- data.frame(ID= 1:(nid*nrep), SCEN = rep(1:nrep,nid), REP = rep(1:nid,each=nrep))
set.seed(1) ; eta1 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var1))
set.seed(2) ; eta2 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var2))
set.seed(3) ; eta3 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var3))
set.seed(4) ; eta4 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var4))
set.seed(5) ; eta5 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var5))
set.seed(6) ; eta6 <- rnorm(n=nrow(samp_shell),mean=0,sd=sqrt(var6))
# samp_ETA <-  mutate(samp_shell,ETA1=eta1,ETA2=eta2,ETA3=eta3,ETA4=eta4,ETA5=eta5,ETA6=eta6)
samp_ETA <-  mutate(samp_shell,ETA1=0,ETA2=0,ETA3=0,ETA4=0,ETA5=0,ETA6=0)

set.seed(10) ; adadraw <- runif(nrow(samp_ETA))
samp <- dplyr::left_join(samp_ETA,s_df,by="SCEN") %>%
  # mutate(ADADRAW=adadraw,ADA=ifelse(CMPD==1&SPEC==2&ADADRAW<prob_monkey_ada,1,ifelse(CMPD==1&SPEC%in%c(3,4)&ADADRAW<prob_rodent_ada,1,0)))
  mutate(ADADRAW=1,ADA=ifelse(CMPD==1&SPEC==2&ADADRAW<prob_monkey_ada,1,ifelse(CMPD==1&SPEC%in%c(3,4)&ADADRAW<prob_rodent_ada,1,0)))

# Testing via # x <- 1 # unlist(samp[x,])
out <- lapply(1:nrow(samp), function(x){
  print(x)
  parm  <- par_func(TVCL=samp[x,]$TVCL,TVV1=samp[x,]$TVV1,TVQ=samp[x,]$TVQ,TVV2=samp[x,]$TVV2,TVVMAX=samp[x,]$TVVMAX,TVKM=samp[x,]$TVKM,TVKA=samp[x,]$TVKA,
                    ETA1=samp[x,]$ETA1,ETA2=samp[x,]$ETA2,ETA3=samp[x,]$ETA3,ETA4=samp[x,]$ETA4,ETA5=samp[x,]$ETA5,ETA6=samp[x,]$ETA6,
                    KP=samp[x,]$KP,TVQB=samp[x,]$TVQB,TVVB=samp[x,]$TVVB,
                    TV_SAF_KOUT=samp[x,]$TV_SAF_KOUT,TV_SAF_KIN=samp[x,]$TV_SAF_KIN,TV_EFF_KOUT=samp[x,]$TV_EFF_KOUT,TV_EFF_KIN=samp[x,]$TV_EFF_KIN,
                    SAF_MAX=samp[x,]$SAF_MAX,SAF_C50=samp[x,]$SAF_C50,SAF_GAM=samp[x,]$SAF_GAM,
                    EFF_MAX=samp[x,]$EFF_MAX,EFF_C50=samp[x,]$EFF_C50,EFF_GAM=samp[x,]$EFF_GAM,
                    SPEC=samp[x,]$SPEC,
                    BW=samp[x,]$BW,
                    ADA=samp[x,]$ADA,
                    ADASTART=ada_onset) # If ADA=1, ADA effect kicks in onset of time, currently no IIV
  inits <- init_func(parm)
  if(unlist(samp[x,])[["ROUTE"]]==1) evnt  <- dose_func(cmt=1,
                                                        value=samp[x,]$DOSE_MG,
                                                        # tinf=2, # careful for infusion dose in cmt=3, for bolus dose in cmt=1
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  if(unlist(samp[x,])[["ROUTE"]]!=1) evnt  <- dose_func(cmt=5,
                                                        value=samp[x,]$DOSE_MG * samp[x,]$TVF,
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  
  # times  <- sort(c(1,seq(0,29*7*24,1)))
  times  <- c(0,1*7*24,27*7*24,28*7*24)
  
  data.frame(deSolve::lsoda(inits,times,des_func,parm,events=list(data=evnt)),
             ID=unlist(samp[x,])[["ID"]]
  )
})
out <- do.call(rbind,out)

outsel_hed <- dplyr::select(out,time,C1,AUC=A4,CB,SAFETY=A7,EFFICACY=A8,ID) %>% 
  dplyr::left_join(dplyr::select(samp,ID:REP,SPEC:CMPD,BW,DOSE_MG,ADA),by="ID") %>%
  mutate(SPEC=ifelse(SPEC==1,"Human",ifelse(SPEC==2,"Monkey",ifelse(SPEC==3,"Rat",ifelse(SPEC==4,"Mouse","ERROR")))),
         CMPD=as.factor(ifelse(CMPD==1,"mAb","SM")))

outsel_hed_ss_start <- subset(outsel_hed,NDOS!=1&time==27*7*24)
outsel_hed_ss_end <- subset(outsel_hed,NDOS!=1&time==28*7*24) %>% dplyr::select(ID,AUC_END=AUC)
outsel_hed_ss <- dplyr::left_join(outsel_hed_ss_start,outsel_hed_ss_end,by="ID") %>% mutate(AUCSS=AUC_END-AUC)

pl_hed_sd <- ggplot(data=subset(outsel_hed,NDOS==1&time==28*24*7&DOSE_MGKG<15))+
    geom_line(aes(x=DOSE_MGKG,y=AUC))+
    geom_hline(yintercept = auc_noael_sd/10,lty=2) + 
    geom_vline(xintercept = 10.9,lty=2) +
    scale_x_continuous(breaks=seq(0,100,1))+
    labs(x="Dose (mg/kg)",y="AUC (h*?g/mL)",subtitle="mAb: Human vs NOAEL AUCinf,sd")+
    theme_lapp()+theme(legend.position="bottom")

pl_hed_w1 <- ggplot(data=subset(outsel_hed,NDOS!=1&time==1*24*7&DOSE_MGKG<15))+
  geom_line(aes(x=DOSE_MGKG,y=AUC))+
  geom_hline(yintercept = auc_noael_w1/10,lty=2) + 
  geom_vline(xintercept = 8.8,lty=2) +
  scale_x_continuous(breaks=seq(0,100,1))+
  labs(x="Dose (mg/kg)",y="AUC (h*?g/mL)",subtitle="mAb: Human vs NOAEL AUC0-7d")+
  theme_lapp()+theme(legend.position="bottom")

pl_hed_ss <- ggplot(data=subset(outsel_hed_ss,DOSE_MGKG<15))+
  geom_line(aes(x=DOSE_MGKG,y=AUCSS))+
  geom_hline(yintercept = auc_noael_ss/10,lty=2) + 
  geom_vline(xintercept = 7.2,lty=2) +
  scale_x_continuous(breaks=seq(0,100,1))+
  labs(x="Dose (mg/kg)",y="AUC (h*?g/mL)",subtitle="mAb: Human vs NOAEL AUCss,qw")+
  theme_lapp()+theme(legend.position="bottom")

pl_mabel_w1 <- ggplot(data=subset(outsel_hed,NDOS!=1&time==1*24*7&DOSE_MGKG<3))+
  geom_line(aes(x=DOSE_MGKG,y=AUC/168))+
  geom_hline(yintercept = 1,lty=2) +
  geom_hline(yintercept = 2,col="dodgerblue2",lty=2) +
  geom_vline(xintercept = 0.8,lty=2) +
  geom_vline(xintercept = 1.1,col="dodgerblue2",lty=2) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  # geom_vline(xintercept = 8.8,lty=2) +
  scale_x_continuous(breaks=seq(0,5,0.2))+
  labs(x="Dose (mg/kg)",y="Caverage (?g/mL)",subtitle="mAb: MABEL - First week average")+
  theme_lapp()+theme(legend.position="bottom")

pl_mabel_ss <- ggplot(data=subset(outsel_hed_ss,DOSE_MGKG<3))+
  geom_line(aes(x=DOSE_MGKG,y=AUCSS/168))+
  geom_hline(yintercept = 1,lty=2) +
  geom_hline(yintercept = 2,col="dodgerblue2",lty=2) +
  geom_vline(xintercept = 0.8,lty=2) +
  geom_vline(xintercept = 1.1,col="dodgerblue2",lty=2) +
  scale_y_log10()+  annotation_logticks(sides = "l") +
  # geom_vline(xintercept = 8.8,lty=2) +
  scale_x_continuous(breaks=seq(0,5,0.2))+
  labs(x="Dose (mg/kg)",y="Caverage (?g/mL)",subtitle="mAb: MABEL - Steady state week average")+
  theme_lapp()+theme(legend.position="bottom")


# Extra - Trainingsmodule NM - Create NM dataset --------------------------

head(nmds_sd_pk)

nmds_sd_pk <- dplyr::select(trans_sd_err,ID,TIME=time,DV=C1_Y,FLAG=C1_BQL,DOSE_MGKG,ROUTE,CMPD,BW,ADA) %>%
  mutate(CMPD=ifelse(CMPD=="mAb",1,2),CMT=2,EVID=0,LDV=ifelse(DV==0,log(0.05)/2,log(DV)),AMT=0)

nmds_sd_drec <- subset(nmds_sd_pk,!duplicated(ID)) %>% 
  mutate(TIME=0,DV=0,FLAG=0,CMT=ifelse(ROUTE==1,2,1),EVID=1,LDV=0,AMT=DOSE_MGKG*BW)

nmds_sd <- rbind(nmds_sd_pk,nmds_sd_drec) %>% 
  dplyr::select(ID,EVID,AMT,TIME,DV,LDV,CMT,FLAG,DOSE_MGKG:ADA) %>%
  mutate(SPEC=ifelse(BW==3,2,ifelse(BW==0.3,3,ifelse(BW==0.02,4,1))))

nmds_sd <- nmds_sd[order(nmds_sd$ID,nmds_sd$TIME), ]
nmds_sd$RECN <- 1:nrow(nmds_sd)

write.csv(nmds_sd,"../../EXPL.ANALYSIS/SDPK/nmds_sd.csv", row.names = FALSE, quote = FALSE, na = ".")


# NMDS PKPD ---------------------------------------------------------------

nmds_pkpd_mut <- read.csv("trans_pkpd_q2d_long.csv") %>%
  mutate(CMT=ifelse(TYPE=="Cplasma",2,
                    ifelse(TYPE=="Safety_Biomarker",3,
                           ifelse(TYPE=="Cbrain",4,5))),
                    EVID=0,
                    LDV=ifelse(DV==0,log(0.05)/2,log(DV)),
                    LDV2=ifelse(DV==0,log(0.05)/2,ifelse(TYPE%in%c("Cplasma","Cbrain"),log(DV),DV)),AMT=0,
         TYPE=ifelse(TYPE=="Cplasma",2,
                    ifelse(TYPE=="Safety_Biomarker",3,
                           ifelse(TYPE=="Cbrain",4,5))),
         FLAG=ifelse(DV==0,1,0))

nmds_pkpd_drec <- subset(nmds_pkpd_mut,!duplicated(ID)) %>% 
  mutate(TIME=0,DV=0,FLAG=0,TYPE=1,CMT=ifelse(ROUTE==1,2,1),EVID=1,LDV=0,LDV2=0,AMT=DOSE_MGKG*BW)

nmds_pkpd_safbsl <- subset(nmds_pkpd_mut,CMT==3) %>% subset(!duplicated(ID)) %>% dplyr::select(ID,SAFBSL = BSL)
nmds_pkpd_effbsl <- subset(nmds_pkpd_mut,CMT==5) %>% subset(!duplicated(ID)) %>% dplyr::select(ID,EFFBSL = BSL)

nmds_pkpd <- rbind(nmds_pkpd_mut,nmds_pkpd_drec) %>% 
  dplyr::select(ID,EVID,AMT,TIME,DV,LDV,LDV2,CMT,TYPE,FLAG,DOSE_MGKG:ADA) %>%
  mutate(SPEC=ifelse(BW==3,2,ifelse(BW==0.3,3,ifelse(BW==0.02,4,1)))) %>%
  dplyr::left_join(nmds_pkpd_safbsl) %>%
  dplyr::left_join(nmds_pkpd_effbsl)

nmds_pkpd <- nmds_pkpd[order(nmds_pkpd$ID,nmds_pkpd$TIME,nmds_pkpd$CMT), ]
nmds_pkpd$RECN <- 1:nrow(nmds_pkpd)

write.csv(nmds_pkpd,"../../EXPL.ANALYSIS/PKPD/nmds_pkpd.csv", row.names = FALSE, quote = FALSE, na = ".")


# Learning topic list --------------------------------------------------------------

# Using the data exploration tool, following can be learned

# Key learning point Graham mentioned: variability

# In rat (behalve leren allometrie)
# mab niet linear: punt 3), 6)
# SM is de lineare, + flip flop

# NU gaan we rd -> rat
# predictie naar apen
# wij tonen dan mens sprong (samen)

# 1) Linear PK: dose linearity (NCA DN AUC plot)                        -nca_sd_trans_lin.csv
# 2) Linear PK: allometric scaling (NCA parameter plots)                -nca_sd_trans_lin.csv
# 3) Non linear PK: identify this issue through exploration             -transl_sd_err.csv ; nca_sd_trans.csv
# 4) Extravascular PK: identify flip flop, bioavailability              -trans_lin_sd_po_err.csv ; nca_sd_trans_lin_po.csv
# 5) Extravascular PK: time points matter - missed Cmax                 -trans_lin_sd_po_err.csv ; nca_sd_trans_lin_po.csv
# 6) ADA mAbs: identify                                                 -transl_sd_err.csv ; nca_sd_trans.csv
# 7) Step by step development: PKPD in 2 compounds in rat and monkey    -trans_pkpd_q2d_long.csv
# 8) At end, show true model input: IPRED vs DV ->                      -pkpd_q2d_trans_n10.csv
#     model scheme + parameters
#     e.g. no variability in drug effect and still difficult to pick up due to error + IIV
#     currently no species differences, this can be squeezed in to show something specifically but rapidly complex


# To Do -------------------------------------------------------------------

# Poor cross react in rat : draai aan EC50 - mAb doet het niet in rat
# Affinity: human=monkey<<<rat voor de mAb

# Maak dan ook een typical sim om te tonen hoe het model werkt (en dan naar mens + IIV band)

# Dit kan meteen in een interne case study gegoten worden 
# (wel KP verhaal naar iets anders omdat dit lijkt op preclin DD course)
# Kunnen we mensen exploratie app laten gebruiken zonder te mogen fitten -> dat is kern van die case dan
# Idee: wat als de groepjes me dan 1x mogen een study design kunnen geven ->
# Ik simuleer die on the fly en geef ze de output na een uur of zo
# Zo kunnen ze verder aan de slag om beter te begrijpen :)
# bvb ze krijgen rat -> ze geven me een monkey design -> daar krijgen ze results van met hun gekozen dose/tps..
# dan moeten ze met rat + aap de translatie naar mens maken
