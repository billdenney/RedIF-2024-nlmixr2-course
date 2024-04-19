#
# nlmxr model code for a 
# two compartment model with first order absorption and
# and effect compartment
# 
# f fixed to 0.6 --> not sure how to implement in closed form
#
# allometric scaling on all clearance and volume parameters
# 
#
# prop error on cp
#
#################################################################################

cmt2_ode <- function() {
  ini({
    
    # use mu referencing 
    tvka     <- log(1)
    tvcl     <- log(0.4)   # log Cl (L/h)
    tvvcent  <- log(3)     # log Vcent (L)
    tvQ      <- log(0.8)     # log Q (L/h)
    tvvperip  <- log(9)     # log Vcent (L)
    tvfdepot <-  fix(0.6)     # bioavailability
    
    #exponents for allometric scaling
    ascl     <- fix(0.75)    # exponent all. scaling of CL and Q
    asv      <- fix(1)       # exponent all. scaling of CL and Q
    
    # iiv 
    eta.ka     ~ fix(0) # eta on ka
    eta.cl     ~ 0.1 # eta on cl as variance
    eta.vcent  ~ 0.04 # eta on Vcent as variance
    eta.q      ~  fix(0) # eta on q
    eta.perip  ~  fix(0) # eta on Vperh
    
    # residual error
    eps.pkprop <-  0.3
  })
  model({

    
    # PK parameter
    ka    = exp(tvka + bw_scl + eta.ka)
    CL    = exp(tvcl + bw_scl*ascl + eta.cl)
    Vcent = exp(tvvcent + bw_scl*asv + eta.vcent)
    Q     = exp(tvQ + bw_scl*ascl + eta.q)
    Vcper = exp(tvvperip  + bw_scl*asv + eta.perip)
    Fdepot = tvfdepot 
    
    kcp   = Q/Vcent
    kpc   = Q/Vcper
    kel   = CL/Vcent
    
    
 
    d/dt(depot)    =  -ka * depot
    f(depot) = Fdepot
    d/dt(central)  =  ka  * depot - (kcp+kel) * central + kpc * peripher
    d/dt(peripher) =  kcp * central - kpc * peripher

    
    CP = central/Vcent

    
    # additive error to start with
    CP ~ prop(eps.pkprop)
  })
}

### closed form --> dont know how to use F in this case 

cmt2 <- function() {
  ini({
    
    # use mu referencing 
    tvka     <- log(0.001)
    tvcl     <- log(0.4)   # log Cl (L/h)
    tvvcent  <- log(3)     # log Vcent (L)
    tvQ      <- log(0.8)     # log Q (L/h)
    tvvperip  <- log(9)     # log Vcent (L)
    tvfdepot <-  fix(0.6)     # bioavailability
    
    #exponents for allometric scaling
    ascl     <- fix(0.75)    # exponent all. scaling of CL and Q
    asv      <- fix(1)       # exponent all. scaling of CL and Q
    
    # iiv 
    eta.ka     ~ fix(0) # eta on ka
    eta.cl     ~ 0.1 # eta on cl as variance
    eta.vcent  ~ 0.04 # eta on Vcent as variance
    eta.q      ~  fix(0) # eta on q
    eta.perip  ~  fix(0) # eta on Vperh
    
    # residual error
    eps.pkprop <-  0.3
  })
  model({
    
    
    # PK parameter
    KA    = exp(tvka + bw_scl + eta.ka)
    CL    = exp(tvcl + bw_scl*ascl + eta.cl)
    V2 = exp(tvvcent + bw_scl*asv + eta.vcent)
    Q     = exp(tvQ + bw_scl*ascl + eta.q)
    V3 = exp(tvvperip  + bw_scl*asv + eta.perip)
    Fdepot = tvfdepot 



    linCmt(KA, CL, V2, Q, V3) ~ prop(eps.pkprop) 


})
}