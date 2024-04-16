##### 
#  Title: Simulate human PK/PD dataset for fictive small molecule.
# 
#     Model was also used to simulate pre-clinical data. 
#     
#     
#
#
#
# Libraries ---------------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(pmxsimtools)
library(dplyr)

home <- getwd()


dir.create("datasets")


# Redif project ----------------------------------------------------------------

# set and  create folders ------------------------------------------------------


mod_fold <- "models_r/"


# ---- PKPD simulations settings  ----------------------------------------------
# 
#  Simulate everything in the original setting but just provide the small molecule data only.
#


# fictive PD

# Safety driven by Cplasma -->  example albumin drop --> not in animals now in humans
# Efficacy driven by Cplasma -->  example albumin drop --> not in animals now in humans


# PKPD - Simulate  clinical ------------------------------
#
#
#
#
# User input
source(paste0(mod_fold,"ode_translation_PKPD_red.r"))


# Fictive Safety and Efficacy Turnover
# Note: I call kout TV, because in the model I'll scale with body weight (-0.25), not for Kin (choice at this point)
# Note: I call kin TV, because in the model we'll add IIV

# IIV
varCL <- 0.1     # variance CL
varV1 <- 0.04    # variance V1
varKSA <- 0.1     # variance KIN SAFETY
varKEF <- 0.1     # variance KIN EFFICACY

# Error SD
sd_cplasma_prp_err <- 0.3
sd_cbrain_prp_err  <- 0.4
sd_eff_prp_err     <- 0.15
sd_saf_prp_err     <- 0.1

# LLOQ 
lloq_cplasma <- 0.05
lloq_cbrain  <- 0.05
lloq_eff     <- 0.05
lloq_saf     <- 0.05



# Step 5: human Phase I HV and Phase 2 in patients ------------------------------
#
#  HV comparable to scaled model with increase Ka --> no absorption controlled PK
#
#  Patients:  30 lower CL





p_sm  <- data.frame(TVCL=0.4,
           TVV1=3,
           TVKA=1,
           TVF=0.6,
           TVQ=0.8,
           TVV2=9,
           TVVMAX=0,
           TVKM=1,
           KP = 1    ,
           TVQB=0.02, 
           TVVB = 0.03,
           SAF_MAX = -0.25, 
           SAF_C50 = 10, 
           SAF_GAM = 2, 
           EFF_MAX = -1, 
           EFF_C50 = 2, 
           EFF_GAM = 1,
           TV_SAF_KOUT = 0.0015, 
           SAF_BSL =  40 ,
           TV_EFF_KOUT =0.15 ,
           EFF_BSL =  10 ,
           FAC_PAT_CL  =0.7,
           FAC_PAT_SAF = 0.7,
           FAC_PAT_EFF = 2
           )





nid <- 5

doses <- c(50,100, 150, 200, 300)
# User input HV times
# dosing starts 1 hour after baseline
times_hv <- c(0,c(0,0.25,0.5,0.75,1,1.5,2,3,4,6,8,12,(c(seq(1,13),15,17)*24)-15/60) +1)
times_pat <- c(0, c(0,1,4,(seq(1,31)*24)-15/60) +1)

 
s_df <- expand.grid(SPEC = c(1,5), # HV 1  and Patients 5
                    DOSE_MG = doses, 
                    ROUTE = c(2), # po
                    CMPD = c(2),
                    II = 24,# sm
                    runin = 1
) %>% bind_cols(.,p_sm  ) %>% mutate(SCEN = 1:n()) %>% group_by(SCEN) %>% left_join(.,expand(.,REP =1:nid )) %>% ungroup()
  
samp <- s_df %>% mutate( ETA1 = rnorm(n=nrow(s_df),mean=0,sd=sqrt(varCL)), #CL
                         ETA2 = 0,                                         #ka              
                     
                         ETA3 = rnorm(n=nrow(s_df),mean=0,sd=sqrt(varKSA)),
                         ETA4 = rnorm(n=nrow(s_df),mean=0,sd=sqrt(varKEF)),
                         BW = rnorm(n=nrow(s_df), mean = 78.3, sd = 10.4 ),
                         ID = 1:n()) %>% 
                 group_by(SPEC, DOSE_MG) %>% 
                 mutate(ID2 = sample(c(1:nid), size =  nid, prob = rep(1/nid,nid))) %>% ungroup() %>%
                 mutate(TVCL = case_when(SPEC == 5 ~ TVCL *FAC_PAT_CL, TRUE ~ TVCL),
                        EFF_BSL = case_when(SPEC == 5 ~ EFF_BSL *FAC_PAT_EFF, TRUE ~ EFF_BSL),
                        SAF_BSL = case_when(SPEC == 5 ~ SAF_BSL *FAC_PAT_SAF, TRUE ~ SAF_BSL),
                        NDOS =  case_when(SPEC == 5 ~ 28, TRUE ~ 10 ),
                        TV_SAF_KIN =  TV_SAF_KOUT * SAF_BSL,
                        TV_EFF_KIN =  TV_EFF_KOUT * EFF_BSL)



# Testing via # x <- 1 # unlist(samp[x,])
out <- lapply(1:nrow(samp), function(x){ #x = 1
  print(x)
  parm  <- par_func(TVCL=samp[x,]$TVCL,TVV1=samp[x,]$TVV1,TVQ=samp[x,]$TVQ,TVV2=samp[x,]$TVV2,
                    TVKA=samp[x,]$TVKA,TVVMAX=samp[x,]$TVVMAX,TVKM=samp[x,]$TVKM,
                    KP=samp[x,]$KP,TVQB=samp[x,]$TVQB,TVVB=samp[x,]$TVVB,
                    ETA1=samp[x,]$ETA1,ETA2=samp[x,]$ETA2,ETA3=samp[x,]$ETA3,ETA4=samp[x,]$ETA4,
                    TV_SAF_KOUT=samp[x,]$TV_SAF_KOUT,TV_SAF_KIN=samp[x,]$TV_SAF_KIN,
                    TV_EFF_KOUT=samp[x,]$TV_EFF_KOUT,TV_EFF_KIN=samp[x,]$TV_EFF_KIN,
                    SAF_MAX=samp[x,]$SAF_MAX,SAF_C50=samp[x,]$SAF_C50,SAF_GAM=samp[x,]$SAF_GAM,
                    EFF_MAX=samp[x,]$EFF_MAX,EFF_C50=samp[x,]$EFF_C50,EFF_GAM=samp[x,]$EFF_GAM,
                    SPEC=samp[x,]$SPEC,
                    BW=samp[x,]$BW) # If ADA=1, ADA effect kicks in onset of time, currently no IIV
  inits <- init_func(parm)
  
  if(samp$ROUTE[x]==1) evnt  <- dose_func(cmt=1,
                                                        value=samp[x,]$DOSE_MG,
                                                        # tinf=2, # careful for infusion dose in cmt=3, for bolus dose in cmt=1
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  if(samp$ROUTE[x]!=1) evnt  <- dose_func(cmt=5,
                                                        value=samp[x,]$DOSE_MG * samp[x,]$TVF,
                                                        tau=samp[x,]$II,
                                                        ndose=samp[x,]$NDOS)
  
  
  evnt$time <- evnt$time + samp$runin[x]
### select time for patients and for hv
  times <- if(samp$SPEC[x] == 1){
    
        times_hv
    
  }else{
    
    times_pat
    
  }
  dd <- data.frame(deSolve::lsoda(inits,times,des_func,parm,events=list(data=evnt)),
             ID=unlist(samp[x,])[["ID"]],
             ID2 =  samp$ID2[x],
             ICL =  parm["CL"],
             IV1 =  parm["V1"],
             IV2 =  parm["V2"],
             IQ =  parm["Q"],
             IKA =  parm["KA"]
            
             
  )
  dd$EVID <- ifelse(dd$time %in% evnt$time ,1,0)
  dd$AMT  <- ifelse(dd$time %in% evnt$time ,samp[x,]$DOSE_MG,0)
return(dd)
  })
out <- do.call(rbind,out)




outsel_human <- dplyr::select(subset(out, AMT <= 0) ,time,EVID,AMT, CP_IPRED = C1, SAF_IPRED=A7,EFF_IPRED=A8,ID, ID2,ICL,IV1, IV2, IQ, IKA) %>%   mutate(
         CP_ERR =  rnorm(n = n(), mean = 0, sd = sd_cplasma_prp_err),
         
         EFF_ERR = rnorm(n = n(), mean = 0, sd = sd_eff_prp_err),
         SAF_ERR = rnorm(n = n(), mean = 0, sd = sd_saf_prp_err)
         ) %>% 
  pivot_longer(cols  = !c(ID, ID2,time,EVID,ICL,IV1, IV2, IQ, IKA, AMT), 
                            names_sep = "_",
                            names_to = c("DVID", ".value")) %>% mutate(DV =  exp(log(IPRED) + ERR), 
                                                                       BLOQ = case_when(DV <  0.05 ~ 1,
                                                                                                          TRUE~  0),
                                                                       EVID = case_when(DV <  0.05 ~ 1,
                                                                                                 TRUE~  0),
                                                                       DV = case_when(BLOQ == 1 ~ 0,
                                                                                        TRUE~   DV)
                                                                        )
outsel_human <- outsel_human  %>% left_join(. , select(samp,ID,SPEC, ROUTE, CMPD,BW,DOSE_MG),by="ID")

outsel_human_dose <- out %>% filter( AMT > 0)  %>% select(ID, ID2,time,EVID,ICL,IV1, IV2, IQ, IKA, AMT) %>%  
                             mutate(DVID = "DEPOT", DV = 0, IPRED = 0, ERR = 0, BLOQ = 0)  %>% 
                          left_join(. , select(samp,ID,SPEC, ROUTE, CMPD,BW,DOSE_MG),by="ID")



nml_dd <- bind_rows(outsel_human, outsel_human_dose) %>% mutate( PAT = case_when(SPEC == 5 ~ "Patient",
                                                                                      TRUE ~ "HV"),
                                                                 TIME = time, 
                                                                 CMPD = "small molecule",
                                                                 ROUTE = "oral", 
                                                                 time = NULL) %>% 
                                                        select(ID, ID2, TIME, DV, DVID, EVID,BLOQ, 
                                                               AMT, ROUTE, ICL ,  IV1,   IV2 ,   IQ,   IKA ,
                                                               BW ,DOSE_MG, PAT) %>%
                                                        arrange(ID, TIME, DVID)
  

### write datasets as rds and csv

#HV 50 150 300
write.csv(subset(nml_dd, PAT == "HV" & DVID %in%  c("CP","DEPOT")  %in% c(50,150,300), -c(ID2,ICL ,  IV1,   IV2 ,   IQ,   IKA )),"../datasets/nmds_hv_pk_3doses.csv", row.names = FALSE, quote = FALSE, na = ".")
write.csv(subset(nml_dd, PAT != "HV" & DVID %in%  c("CP","DEPOT") %in% c(150,300), -c(ID2,ICL ,  IV1,   IV2 ,   IQ,   IKA )),"../datasets/nmds_patient_pk.csv", row.names = FALSE, quote = FALSE, na = ".")

# Patients only two doses
write.csv(subset(nml_dd, PAT == "HV" & !DVID %in%  c("CP") & DOSE_MG %in%c(50,150,300), -c(ID2 )),"../datasets/nmds_hv_pkpd_3doses.csv", row.names = FALSE, quote = FALSE, na = ".")
write.csv(subset(nml_dd, PAT != "HV" & !DVID %in%  c("CP") & DOSE_MG %in% c(150,300), -c(ID2)),"../datasets/nmds_patient_pkpd.csv", row.names = FALSE, quote = FALSE, na = ".")


