# Function to define parameters
par_func <- function(TVCL,TVV1,TVQ,TVV2,TVVMAX,TVKM,ETA1,ETA2,ETA3,ETA4,BW,SPEC,TVKA,KP,TVQB,TVVB,
                     TV_SAF_KIN,TV_SAF_KOUT,TV_EFF_KIN,TV_EFF_KOUT,
                     SAF_MAX,SAF_C50,SAF_GAM,EFF_MAX,EFF_C50,EFF_GAM){ 
  
  EFF_KOUT <- TV_EFF_KOUT *((BW/70)**(-0.25))
  SAF_KOUT <- TV_SAF_KOUT *((BW/70)**(-0.25))
  
  SAF_KIN     <-  TV_SAF_KIN * exp(ETA3)
  EFF_KIN     <-  TV_EFF_KIN * exp(ETA4)
  
  KA     <-  TVKA 
  
  CL     <-  TVCL *((BW/70)**0.75) * exp(ETA1) 
  V1     <-  TVV1 *(BW/70) * exp(ETA2) 
  Q      <-  TVQ *((BW/70)**0.75)     
  V2     <-  TVV2 *(BW/70)
  VMAX   <-  TVVMAX *((BW/70)**0.75)
  
  QB     <-  TVQB *((BW/70)**0.75)     
  VB     <-  TVVB *(BW/70)
  
  KM     <-  TVKM
  
  K10    <-  CL/V1
  K12    <-  Q/V1		
  K21    <-  Q/V2
  
  K16  <-  KP*QB/V1
  K61  <-  QB/VB
  

  unlist(mget(ls()))
  
}

# Function to define model
des_func <- function(t,y,p){
  with(as.list(c(p,y)), {
    C1 <- A1/V1
    CB <- A6/VB
   
    DADT1 <-  KA*A5 + A3 -(K10+K12+K16)*A1 + K21*A2 + K61*A6 - VMAX *C1/(KM + C1) 
    DADT2 <-        K12*A1   -K21*A2
    DADT3 <- 0
    DADT4 <- C1
    DADT5 <- -KA*A5
    DADT6 <- K16*A1   -K61*A6
    DADT7 <- SAF_KIN*(1+SAF_MAX*(C1**SAF_GAM)/((SAF_C50**SAF_GAM)+(C1**SAF_GAM))) - SAF_KOUT*A7
    DADT8 <- EFF_KIN*(1+EFF_MAX*(CB**EFF_GAM)/((EFF_C50**EFF_GAM)+(CB**EFF_GAM))) - EFF_KOUT*A8
    list(c(DADT1, DADT2, DADT3, DADT4, DADT5, DADT6, DADT7, DADT8),c(C1=C1,CB=CB))
  })
}

# Function to define initials
init_func <- function(pars){
  init <- c(
    A1 = 0,
    A2 = 0,
    A3 = 0,
    A4 = 0,
    A5 = 0,
    A6 = 0,
    A7 = pars['SAF_KIN']/pars['SAF_KOUT'],
    A8 = pars['EFF_KIN']/pars['EFF_KOUT']
  )
  names(init) <- paste0("A",1:length(init))
  return(init)
}

dose_func<- function (cmt, value, tinf, tau, ndose, times) 
{
  if (missing(times)) {
    timing <- seq(0, (ndose - 1) * tau, tau)
  }
  else {
    timing <- times
  }
  dose <- data.frame(var = paste0("A", cmt), time = timing, 
                     value = value, method = "add")
  if (!missing(tinf)) {
    dose2 <- data.frame(var = paste0("A", cmt), time = timing + 
                          tinf, value = 0, method = "rep")
    dose <- rbind(dose, dose2)
    dose <- dose[order(dose$time), ]
    dose$value[dose$value != 0] <- dose$value[dose$value != 
                                                0]/tinf
  }
  return(dose)
}
