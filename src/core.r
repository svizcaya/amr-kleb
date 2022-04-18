## Core file includes all functions: fits, mcmc sampling, mcmc solutions, sensitivity, plots

require("deSolve")
require("ggplot2")
require("reshape2")
require("gridExtra")
require("FME")
require("Rcgmin")

parallel.fits.f = function(nu.inf.sens.parallel.value,       
                           nu.inf.sensitivity.logical = FALSE,
                           nu.inf.sens.parallel,
                           epsilon.equivalent.sensitivity, 
                           init.solution,
                           quiere.sensitivity,    
                           fit.incl.epsilon = FALSE,
                           fit.excl.epsilon = FALSE,
                           baseline.epsilon,
                           mcmc.compute,
                           mcmc.solve,
                           mcmc.plot,
                           optimize.epsilon,
                           fitfolder,
                           amc.scenarios,
                           Figure = "NADA"){
    graphics.off()

    ## -Estimated(from data) risk in a fully susceptible population exponentially distributed time to colonization  
    r1.last.data   <<- 802/21147
    t.r1.last.data <<- 2017.16-2000

    ## steer flags to control optimization execution
    plot.fits = FALSE
    quiere.fits = FALSE
    
    if(fit.incl.epsilon == TRUE | fit.excl.epsilon == TRUE ) {quiere.fits = TRUE}
    
    nu.inf.sens = nu.inf.sens.parallel.value
    
    source("handy_functions.r")
    
    tiffname1 = "tiffname1_aux"
    tiffname2 = "tiffname2_aux"
    numbers_file = "numbers_file_aux"

    ## just remainders to expect long runs
    
    if(fit.incl.epsilon == TRUE){message("Setting to run Fit including epsilon MCMC")}
    if(fit.excl.epsilon == TRUE){message("Fit excluding epsilon MCMC")}
    if(mcmc.compute == TRUE){message("compute MCMC excluding epsilon MCMC")}

    ## configure settings to comupute and plot figures
    ## -figures group 3: scenarios
    ## --3A scenarios of overall consumption     
    if(Figure == "3A"){
        message("DRAW FIGURE 3A with MCMC")
        figure.label =  "Figure3A"
        beta.1.h.scenarios = FALSE
        amc.scenarios = TRUE
        amc.sc.onlyrestricted.st = FALSE
        tiffname1 = "Figure3A_MCMC.tiff"
        tiffname2 = "Figure3A_MCMC_axis.tiff"
        numbers_file = paste(fitfolder, "Figure_3A_numbers_MCMC", sep = "/")
        tiffname1 = paste(fitfolder, tiffname1, sep = "/")
        tiffname2 = paste(fitfolder, tiffname2, sep = "/")
        ylimite = c(0,50)
    }
    ## --3B scenarios of restricted type consumption     
    if(Figure == "3B"){
        message("DRAW FIGURE 3B with MCM")
        figure.label = "Figure3B"
        beta.1.h.scenarios = FALSE
        amc.scenarios = TRUE
        amc.sc.onlyrestricted.st = TRUE
        tiffname1 = "Figure3B_MCMC.tiff"
        tiffname2 = "Figure3B_MCMC_axis.tiff"
        tiffname1 = paste(fitfolder, tiffname1, sep = "/")
        tiffname2 = paste(fitfolder, tiffname2, sep = "/")
        numbers_file = paste(fitfolder, "Figure_3B_numbers", sep = "/")
        ylimite = c(0,50)        
    }
    ## --3C scenarios of hospital transmission
    if(Figure =="3C"){
        message("DRAW FIGURE 3C MCMC")
        figure.label = "Figure3C"
        beta.1.h.scenarios = TRUE
        amc.scenarios = FALSE
        amc.sc.onlyrestricted.st = FALSE
        tiffname1 = "Figure3C_MCMC.tiff"
        tiffname2 = "Figure3C_axis.tiff"
        tiffname1 = paste(fitfolder, tiffname1, sep = "/")
        tiffname2 = paste(fitfolder, tiffname2, sep = "/")
        plot.fits = FALSE
        numbers_file = paste(fitfolder, "Figure_3C_numbers", sep = "/")
        ylimite = c(0,50)
    }
    ## -figures group 2: sensitivity on epsilon
    if(Figure == "2"){
        message("DRAW FIGURE 2 MCMC")
        epsilons.share = seq(0, 0.60, by = 0.15)                
        figure.label = "Figure2"
        amc.sc.onlyrestricted.st = FALSE
        beta.1.h.scenarios = FALSE
        amc.scenarios = FALSE
        tiffname1 = "Figure2_MCMC.tiff"
        tiffname2 = "Figure2_axis.tiff"
        tiffname1 = paste(fitfolder, tiffname1, sep = "/")
        tiffname2 = paste(fitfolder, tiffname2, sep = "/") 
        ylimite = c(0,20)
        xlimite = c(2000,2020)
    }


    ## hardcoded input values
    ## -number of individuals (normalized) 
    canton.population.init =  100000
    
    ## -unit of time is years, origin 2000
    initime = 0
    endtime = 17.5
        
    init.hospitalized.fraction = 0.001723956

    ## -prevalence observed data points
    prevalence.resistance.ceftriaxon.dir.inpatient =  data.frame(time = 4.5:17.5,
                                                                 prop.r1.h = c(
                                                                     0.0,
                                                                     1.4,
                                                                     7.9,
                                                                     1.3,
                                                                     2.2,
                                                                     1.3,
                                                                     4.9,
                                                                     3.9,
                                                                     2.5,
                                                                     4.5,
                                                                     6.4,
                                                                     4.6,
                                                                     8.4,
                                                                     10.4)/100,
                                                                 numerador = c(
                                                                     0,
                                                                     4,
                                                                     24,
                                                                     4,
                                                                     7,
                                                                     5,
                                                                     17,
                                                                     20,
                                                                     12,
                                                                     26,
                                                                     38,
                                                                     29,
                                                                     59,
                                                                     10),
                                                                 denominador = c(
                                                                     208,
                                                                     278,
                                                                     305,
                                                                     300,
                                                                     325,
                                                                     374,
                                                                     344,
                                                                     508,
                                                                     488,
                                                                     579,
                                                                     597,
                                                                     627,
                                                                     701,
                                                                     96)
                                                                 )
    
    prevalence.resistance.ceftriaxon.dir.outpatient = data.frame(time = 4.5:17.5,
                                                                 prop.r1.c = c(
                                                                     NA,
                                                                     NA,
                                                                     NA,
                                                                     2.7,
                                                                     1.17,
                                                                     0,
                                                                     1.96,
                                                                     1.60,
                                                                     2.84,
                                                                     5.02,
                                                                     2.23,
                                                                     1.56,
                                                                     2.36,
                                                                     NA)/100,
                                                                 numerador =  c(
                                                                     NA,
                                                                     NA,
                                                                     NA,
                                                                     4.00,
                                                                     2.00,
                                                                     0.00,
                                                                     3.00,
                                                                     3.00,
                                                                     6.00,
                                                                     11.00,
                                                                     6.00,
                                                                     5.00,
                                                                     7.00,
                                                                     NA),
                                                                 denominador =  c(
                                                                     NA,
                                                                     NA,
                                                                     NA,
                                                                     146.00,
                                                                     171.00,
                                                                     182.00,
                                                                     153.00,
                                                                     188.00,
                                                                     211.00,
                                                                     219.00,
                                                                     269.00,
                                                                     321.00,
                                                                     296.00,
                                                                     NA
                                                                 ))
    




    ## -estimate poisson confidence level intervals
    prevalence.resistance.ceftriaxon.dir.inpatient = subset(prevalence.resistance.ceftriaxon.dir.inpatient, numerador !=0)
    prevalence.resistance.ceftriaxon.dir.outpatient = subset(prevalence.resistance.ceftriaxon.dir.outpatient, numerador !=0)

    prevalence.resistance.ceftriaxon.dir.inpatient$prop.r1.h = prevalence.resistance.ceftriaxon.dir.inpatient$numerador/prevalence.resistance.ceftriaxon.dir.inpatient$denominador
    prevalence.resistance.ceftriaxon.dir.inpatient$ci.upper = qchisq( df = 2*(prevalence.resistance.ceftriaxon.dir.inpatient$numerador + 1) , p = 1 - 0.05/2 ) / 2 / prevalence.resistance.ceftriaxon.dir.inpatient$denominador
    prevalence.resistance.ceftriaxon.dir.inpatient$ci.lower = qchisq( df = 2*prevalence.resistance.ceftriaxon.dir.inpatient$numerador , p = 0.05/2 ) / 2 / prevalence.resistance.ceftriaxon.dir.inpatient$denominador
    prevalence.resistance.ceftriaxon.dir.inpatient$place = "Hospitals"    
    
    prevalence.resistance.ceftriaxon.dir.outpatient$prop.r1.c = prevalence.resistance.ceftriaxon.dir.outpatient$numerador/prevalence.resistance.ceftriaxon.dir.outpatient$denominador
    prevalence.resistance.ceftriaxon.dir.outpatient$ci.upper = qchisq( df = 2*(prevalence.resistance.ceftriaxon.dir.outpatient$numerador + 1) , p = 1 - 0.05/2 ) / 2 / prevalence.resistance.ceftriaxon.dir.outpatient$denominador
    prevalence.resistance.ceftriaxon.dir.outpatient$ci.lower = qchisq( df = 2*prevalence.resistance.ceftriaxon.dir.outpatient$numerador , p = 0.05/2 ) / 2 / prevalence.resistance.ceftriaxon.dir.outpatient$denominador
    prevalence.resistance.ceftriaxon.dir.outpatient$place = "Community" 

    prevalence.resistance.ceftriaxon.dir.inpatient$std  = sqrt(prevalence.resistance.ceftriaxon.dir.inpatient$prop.r1.h*(1-prevalence.resistance.ceftriaxon.dir.inpatient$prop.r1.h)/prevalence.resistance.ceftriaxon.dir.inpatient$denominador)    
    prevalence.resistance.ceftriaxon.dir.outpatient$std = sqrt(prevalence.resistance.ceftriaxon.dir.outpatient$prop.r1.c*(1-prevalence.resistance.ceftriaxon.dir.outpatient$prop.r1.c)/prevalence.resistance.ceftriaxon.dir.outpatient$denominador)

    write.table(prevalence.resistance.ceftriaxon.dir.inpatient, "prevalence_resistance_ceftriaxon_dir_inpatient")
    write.table(prevalence.resistance.ceftriaxon.dir.outpatient, "prevalence_resistance_ceftriaxon_dir_outpatient")
        
    ## parameters for the differential equations initializsation
    parameters0 =  list(
        theta.h = NA,
        mu.discharge = 36.5,  #10 days
        mu.death = 0,       
        nu.susc = 3,        
        nu.inf = 2,
        beta.1.h = 12,
        beta.2.h = NA, ## .2 indicates carbapenem resistance. Curretly inactive
        beta.1.c = 0.8,
        beta.2.c = NA,
        epsilon.1.c.slope = epsilon.equivalent.to.slope.f(epsilon.equivalent.sensitivity),
        epsilon.share.slope = 0,
        epsilon.2.c = 0,
        Tau.h.1 = NA,
        Tau.h.2 = NA,
        Tau.h.3 = NA,
        Tau.c.1 = NA,
        Tau.c.2 = NA,
        Tau.c.3 = NA,
        tau.1 = 8/(365.25),#treatment duration
        tau.2 = 5/(365.25),
        tau.3 = 8/(365.25),
        delta = 1/12,## delay - treatment to clearance   ## NOTE: changed during revision 
        omega = 0.5,## relative chance for quasi-neutral antimicrobials
        lambda = 0.80,## proportion who clearn colonization upon treatment        
        alpha.h = - log(1-0.3)/((30.25*3)/365.25), ## 30% in three months
        alpha.c = - log(1-0.3)/((30.25*3)/365.25),
        epsilon.1.h = 0, ## no external force directly into hospitals, i.e. comes through the community
        epsilon.2.h = 0,
        ## init values log transformed equations
        S.normal.init = log( canton.population.init*c(init.hospitalized.fraction, (1-init.hospitalized.fraction)) - c(0,2) ),
        S.amplificado.init       = c(-50,-50),        
        r1.normal.init           = c(-50,-50),
        r1.amplificado.init      = c(-50,-50),
        r1.amplificado.beta.init = c(-50,log (2)),
        r2.normal.init           = c(-50,-50),
        r2.amplificado.init      = c(-50,-50),
        r2.amplificado.beta.init = c(-50,-50),
   
        kleb.prob.h             = 0.07*0.026, ## fraction of antibiotic treatments related to a kleb pneu infection (7% (literature);2.6% invasive from data inpatients)  
        kleb.prob.c             = 0.07*0.019, ## fraction of antibiotic treatments related to a kleb pneu infection (7% (literature);2.6% invasive from data outpatients)  
        phi.h                   = 0,          ## spontaneous clearance in hospital
        phi.c                   = 325.25/30   ## spontaneous clearance of infection
    )

    ## fix nu.inf unelss sensitivity 
    if( nu.inf.sensitivity.logical == TRUE ){parameters0$nu.inf = nu.inf.sens}else{parameters0$nu.inf = 2} 
    epsilon.equivalent.exogeno = epsilon.equivalent.sensitivity
    ## define auxiliary stepwise function
    aux.func.projections =  function(t, base, cut.time = 19, update.value) ifelse(t<= cut.time, 1, 0)*base + ifelse(t<= cut.time,0,1)*as.numeric(update.value)
    ## select parameters to set free for optimization 
    if(optimize.epsilon == TRUE) {parameters.candidates = parameters0[c("beta.1.h", "beta.1.c", "epsilon.1.c.slope")]}else{parameters.candidates = parameters0[c("beta.1.h", "beta.1.c")]}


    ## Function that saturates
    saturating.f = function(t, g, t.interval, t.step, low = 1/2, high = 2, last = FALSE, where = 19, factor = 1/2, till = 20) {
        saturating.f.int = function (tt, g, t.interval, t.step, low, high){
            tp = seq(t.interval[1], t.interval[2], by =  t.step)
            g.tp =     g(tp)
            x.max = max(g.tp)
            x.min = min(g.tp)
            gt          = g(tt)
            aux   = unlist(lapply(gt, function(y) max(x.min*low, min(x.max*high,  y))))
            aux
        }
        out.1 = saturating.f.int(t, g = g, t.interval = t.interval, t.step = t.step, low = low, high = high)
        ## the rate to achieve the factor in 2020:
        if(factor ==0) factor = 1e-16
        factor.rate =  -log(factor)/(till-where)
        df          = ifelse( (t >= where) == 1,1,0 )*(exp(factor.rate*(where-t)))
        out.2 = ifelse( (t < where)  == 1,1,0 )*out.1 + ifelse( (t >= where) == 1,1,0 )*(saturating.f.int(where, g, t.interval, t.step, low, high))*df
        if(last == TRUE) {out =  out.2}else{out =  out.1}
        out
    }
    
    tau.1 =  parameters0$tau.1
    tau.2 =  parameters0$tau.2
    tau.3 =  parameters0$tau.3
    delta =  parameters0$delta

    ## -treatment rate functions (fitted from data)
    ## --hospital
    Tau.1.h.f.base =  function(t){ aux.f = function(tt)  {0.09543    + 0.01071*tt}; (ifelse(t>=6.5,1,0)*aux.f(t) + ifelse(t>=6.5,0,1)*aux.f(6.5)) }
    Tau.2.h.f.base =  function(t){ aux.f = function(tt)  {-0.0006567 + 0.0028152*tt}; (ifelse(t>=6.5,1,0)*aux.f(t) + ifelse(t>=6.5,0,1)*aux.f(6.5)) }
    Tau.3.h.f.base =  function(t){ aux.f = function(tt)  {0.203053   + 0.008512*tt}; (ifelse(t>=6.5,1,0)*aux.f(t) + ifelse(t>=6.5,0,1)*aux.f(6.5)) }
    ## --community
    Tau.1.c.f.base =  function(t){ aux.f = function(tt)  {1.207600   -  0.030742*tt}; (ifelse(t>=13,1,0)*aux.f(t) + ifelse(t>=13,0,1)*aux.f(13))}
    Tau.2.c.f.base =  function(t){ aux.f = function(tt)  {-0.4116815 +  0.0888285*tt - 0.0063209*tt*tt  +  0.0001494*tt*tt*tt}; (ifelse(t>=13,1,0)*aux.f(t) + ifelse(t>=13,0,1)*aux.f(13))}
    Tau.3.c.f.base =  function(t){ aux.f = function(tt)  {2.37326    -  0.05214*tt}; (ifelse(t>=13,1,0)*aux.f(t) + ifelse(t>=13,0,1)*aux.f(13))}
    
    Tau.1.h.f = function(t, factor) 1/(tau.1*365.25/saturating.f(g = Tau.1.h.f.base, t = t, t.interval =  c(7.5,15.5), t.step = 1/12, low = 1/10000, high = 500, last = TRUE, where = 19, factor = factor, till = 20)  + tau.1)
    Tau.2.h.f = function(t, factor) 1/(tau.2*365.25/saturating.f(g = Tau.2.h.f.base, t = t, t.interval =  c(7.5,15.5), t.step = 1/12, low = 1/10000, high = 500, last = TRUE, where = 19, factor = factor, till = 20) + tau.2)
    Tau.3.h.f = function(t, factor) 1/(tau.3*365.25/saturating.f(g = Tau.3.h.f.base, t = t, t.interval =  c(7.5,15.5), t.step = 1/12, low = 1/10000, high = 500, last = TRUE, where = 19, factor = factor, till = 20) + tau.3)
    
    Tau.1.c.f = function(t, factor) 2*1/(tau.1*365.25/saturating.f(g = Tau.1.c.f.base, t = t, t.interval =  c(13,16),    t.step = 1/12, low = 1/10000, high = 500, last = TRUE, where = 19, factor = factor, till = 20) + tau.1)
    Tau.2.c.f = function(t, factor) 2*1/(tau.2*365.25/saturating.f(g = Tau.2.c.f.base, t = t, t.interval =  c(13,16),    t.step = 1/12, low = 1/10000, high = 500, last = TRUE, where = 19, factor = factor, till = 20) + tau.2)
    Tau.3.c.f = function(t, factor) 2*1/(tau.3*365.25/saturating.f(g = Tau.3.c.f.base, t = t, t.interval =  c(13,16),    t.step = 1/12, low = 1/10000, high = 500, last = TRUE, where = 19, factor = factor, till = 20) + tau.3)
    
    tau.convert = function(x, duration){ 1/(duration*365.25/x + duration) }
    
    Tau.1.h.f.aux = Tau.1.h.f
    Tau.2.h.f.aux = Tau.2.h.f
    Tau.3.h.f.aux = Tau.3.h.f
    
    Tau.1.c.f.aux = Tau.1.c.f
    Tau.2.c.f.aux = Tau.2.c.f
    Tau.3.c.f.aux = Tau.3.c.f

    tp =  seq(initime, endtime, by = 1/12)

    ## hospitalisation rate (fit from observed data)
    theta.h.f = function(t) (ifelse(t<=11.5,1,0)*((+0.069-0.063)*1/(1+exp(5*(-t + 7.5))) + 0.063) + ifelse(t<=11.5,0,1)*((+0.069-0.0605)*1/(1+exp(-(-t + 11.5))) + 0.0605))*(1+0.00317)

    ## New people in the canton (microcanonical system)
    theta.c =  function(t) {0+0*t}
    j= (1:2) + 1
    
    bacterial.master <- function(parameters,
                                 update.value.kleb.prob.h =  parameters0["kleb.prob.h"],
                                 update.value.kleb.prob.c =  parameters0["kleb.prob.c"],
                                 cut.time =  100,
                                 amc.factor = 1,
                                 beta.1.h.factor = 1,
                                 amc.sc.onlyrestricted = FALSE){
        
        diffEq =  function(t, state, paramters){ with(as.list(c(state, parameters)),{
            
            kleb.prob.func = function(t) {c(aux.func.projections(t, base = kleb.prob.h, cut.time =  cut.time, update.value = update.value.kleb.prob.h),
                                            aux.func.projections(t, base = kleb.prob.c, cut.time =  cut.time, update.value = update.value.kleb.prob.c))}
            
            epsilon.1.c = ifelse(t<=25, epsilon.1.c.slope*t, epsilon.1.c.slope*25)
            
            beta.2.h = beta.1.h
            beta.2.c = beta.1.c
            beta.1.h = ifelse((t>=19), beta.1.h*(exp(--(t-19)*(log(beta.1.h.factor)/(25-19)))), beta.1.h) ## for scenarios of less transmissions within hospitals
            
            if(amc.sc.onlyrestricted == FALSE){ ## scenarios only for restricted antibiotics
                
                Tau.1.h.f.aux = function(t) Tau.1.h.f(t, factor = amc.factor)
                Tau.2.h.f.aux = function(t) Tau.2.h.f(t, factor = amc.factor)
                Tau.3.h.f.aux = function(t) Tau.3.h.f(t, factor = amc.factor)

                Tau.1.c.f.aux = function(t) Tau.1.c.f(t, factor = amc.factor)
                Tau.2.c.f.aux = function(t) Tau.2.c.f(t, factor = amc.factor)
                Tau.3.c.f.aux = function(t) Tau.3.c.f(t, factor = amc.factor)
                
                Tau.2.delayed.h.f = function(t) 1/(1/saturating.f(g = Tau.2.h.f.base, t = t, t.interval =  c(7.5,15.5), t.step = 1/12, low = 1/5, high = 5, factor = amc.factor) )+ tau.2 + delta
                Tau.2.delayed.c.f = function(t) 1/(1/saturating.f(g = Tau.2.c.f.base, t = t, t.interval =  c(13,16),    t.step = 1/12, low = 1/5, high = 5, factor = amc.factor) )+ tau.2 + delta
            }
            
            if( amc.sc.onlyrestricted == TRUE ){
                Tau.1.h.f.aux = function(t) Tau.1.h.f(t, factor = 1)
                Tau.2.h.f.aux = function(t) Tau.2.h.f(t, factor = amc.factor)
                Tau.3.h.f.aux = function(t) Tau.3.h.f(t, factor = 1)
                
                Tau.1.c.f.aux = function(t) Tau.1.c.f(t, factor = 1)
                Tau.2.c.f.aux = function(t) Tau.2.c.f(t, factor = amc.factor)
                Tau.3.c.f.aux = function(t) Tau.3.c.f(t, factor = 1)
                
                Tau.2.delayed.h.f = function(t) 1/(1/saturating.f(g = Tau.2.h.f.base, t = t, t.interval =  c(7.5,15.5), t.step = 1/12, low = 1/5, high = 5, factor = amc.factor) )+ tau.2 + delta
                Tau.2.delayed.c.f = function(t) 1/(1/saturating.f(g = Tau.2.c.f.base, t = t, t.interval =  c(13,16),    t.step = 1/12, low = 1/5, high = 5, factor = amc.factor) )+ tau.2 + delta
            }
            
            ## --join: originally thougth to support field equations approach
            S.normal.h.c            = c(S.normal.h.c1, S.normal.h.c2)
            S.amplificado.h.c       = c(S.amplificado.h.c1, S.amplificado.h.c2)
            r1.normal.h.c           = c(r1.normal.h.c1, r1.normal.h.c2)
            r1.amplificado.h.c      = c(r1.amplificado.h.c1, r1.amplificado.h.c2)
            r1.amplificado.beta.h.c = c(r1.amplificado.beta.h.c1, r1.amplificado.beta.h.c2)
            r2.normal.h.c           = c(r2.normal.h.c1,      r2.normal.h.c2)
            r2.amplificado.h.c      = c(r2.amplificado.h.c1, r2.amplificado.h.c2)
            r2.amplificado.beta.h.c = c(r2.amplificado.beta.h.c1, r2.amplificado.beta.h.c2)
            ## transform back
            S.normal.h.c.exp            = exp(S.normal.h.c)
            S.amplificado.h.c.exp       = exp(S.amplificado.h.c)
            r1.normal.h.c.exp           = exp(r1.normal.h.c)
            r1.amplificado.h.c.exp      = exp(r1.amplificado.h.c)
            r1.amplificado.beta.h.c.exp = exp(r1.amplificado.beta.h.c)
            r2.normal.h.c.exp           = exp(r2.normal.h.c)
            r2.amplificado.h.c.exp      = exp(r2.amplificado.h.c)
            r2.amplificado.beta.h.c.exp = exp(r2.amplificado.beta.h.c)
                                       
            N = S.normal.h.c.exp + S.amplificado.h.c.exp + r1.normal.h.c.exp + r1.amplificado.h.c.exp + r1.amplificado.beta.h.c.exp + r2.normal.h.c.exp + r2.amplificado.h.c.exp + r2.amplificado.beta.h.c.exp
            
            hg.ratio = 1
            
            Tau.2.f.aux.1 = c(Tau.2.h.f.aux, Tau.2.c.f.aux)
            Tau.2.f.aux = function(t)   unlist(lapply(Tau.2.f.aux.1, function(ff) ff(t)))    -  ifelse(t<=16,1,0)*(1/(1/365.25+tau.2) - c(phi.h, phi.c))*r1.amplificado.h.c.exp/N
            Tau.2.f.S   = function(t)   Tau.2.f.aux(t)
            Tau.2.f.r   = function(t)   hg.ratio*Tau.2.f.S(t)

            kleb.restricted.percent =   100*(1/(1/365.25+tau.2) - c(phi.h, phi.c))*r1.amplificado.h.c.exp/(Tau.2.f.aux(t)*N)
            
            prop.S.normal.h.c.exp            = S.normal.h.c.exp/N
            prop.S.amplificado.h.c.exp       = S.amplificado.h.c.exp/N
            prop.r1.normal.h.c.exp           = r1.normal.h.c.exp/N
            prop.r1.amplificado.h.c.exp      = r1.amplificado.h.c.exp/N
            prop.r1.amplificado.beta.h.c.exp = r1.amplificado.beta.h.c.exp/N
            prop.r2.normal.h.c.exp           = r2.normal.h.c.exp/N
            prop.r2.amplificado.h.c.exp      = r2.amplificado.h.c.exp/N
            prop.r2.amplificado.beta.h.c.exp = r2.amplificado.beta.h.c.exp/N
            
            ## prevalence of colonization   
            prop.r1.h.c = (r1.normal.h.c.exp + r1.amplificado.h.c.exp  + r1.amplificado.beta.h.c.exp + r2.normal.h.c.exp + r2.amplificado.h.c.exp+ r2.amplificado.beta.h.c.exp)/N   

            
            ## differential equations system
            dS.normal.h.c            =   c(1,1)*(1/S.normal.h.c.exp)*(c(theta.h.f(t)*S.normal.h.c.exp[2]                         - mu.discharge*S.normal.h.c.exp[1],              - theta.h.f(t)*S.normal.h.c.exp[2]             + theta.c(t) + mu.discharge*S.normal.h.c.exp[1]  )              - mu.death*S.normal.h.c.exp             - (1/c(N[1],N[2]))*(c(beta.1.h, beta.1.c)*(r1.normal.h.c.exp + nu.inf*r1.amplificado.h.c.exp + nu.inf*r1.amplificado.beta.h.c.exp)                                 + c(beta.2.h, beta.2.c)*(r2.normal.h.c.exp + nu.inf*(r2.amplificado.h.c.exp + r2.amplificado.beta.h.c.exp )))*S.normal.h.c.exp                  - c(Tau.1.h.f.aux(t) + Tau.2.f.S(t)[1]  + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t) + Tau.2.f.S(t)[2] + omega*Tau.3.c.f.aux(t))*S.normal.h.c.exp                         + lambda*c( Tau.2.delayed.h.f(t), Tau.2.delayed.c.f(t) )*(r1.normal.h.c.exp + r1.amplificado.beta.h.c.exp)  + c(alpha.h, alpha.c)*(S.amplificado.h.c.exp + r1.normal.h.c.exp + r2.normal.h.c.exp + r1.amplificado.h.c.exp + r2.amplificado.h.c.exp + r1.amplificado.beta.h.c.exp + r2.amplificado.beta.h.c.exp)   - (c(epsilon.1.h, epsilon.1.c) + c(epsilon.2.h, epsilon.2.c))*S.normal.h.c.exp                 + lambda*(c(1,1)/(1/365.25+tau.2) - c(phi.h, phi.c))*r1.amplificado.h.c.exp)
            
            dS.amplificado.h.c       =   c(1,1)*(1/S.amplificado.h.c.exp)*(c(theta.h.f(t)*S.amplificado.h.c.exp[2]               - mu.discharge*S.amplificado.h.c.exp[1],         - theta.h.f(t)*S.amplificado.h.c.exp[2]        + theta.c(t) + mu.discharge*S.amplificado.h.c.exp[1] )          - mu.death*S.amplificado.h.c.exp        - (1/c(N[1],N[2]))*(c(beta.1.h, beta.1.c)*(r1.normal.h.c.exp + nu.inf*r1.amplificado.h.c.exp + nu.inf*r1.amplificado.beta.h.c.exp)                                 + c(beta.2.h, beta.2.c)*(r2.normal.h.c.exp + nu.inf*(r2.amplificado.h.c.exp + r2.amplificado.beta.h.c.exp )))*nu.susc*S.amplificado.h.c.exp     + c(Tau.1.h.f.aux(t) + Tau.2.f.S(t)[1]  + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t) + Tau.2.f.S(t)[2] + omega*Tau.3.c.f.aux(t))*S.normal.h.c.exp                                                                                                                                     - c(alpha.h, alpha.c)* S.amplificado.h.c.exp                                                                                                                                                          - (c(epsilon.1.h, epsilon.1.c) + c(epsilon.2.h, epsilon.2.c))*nu.susc*S.amplificado.h.c.exp )
            
            dr1.normal.h.c           =   c(1,1)*(1/r1.normal.h.c.exp)*(c(theta.h.f(t)*r1.normal.h.c.exp[2]                       - mu.discharge*r1.normal.h.c.exp[1],             - theta.h.f(t)*r1.normal.h.c.exp[2]            + theta.c(t) + mu.discharge*r1.normal.h.c.exp[1] )              - mu.death*r1.normal.h.c.exp            + (1/c(N[1],N[2]))* c(beta.1.h, beta.1.c)*(r1.normal.h.c.exp + nu.inf*r1.amplificado.h.c.exp + nu.inf*r1.amplificado.beta.h.c.exp)*S.normal.h.c.exp                                                                                                                                                                - c(Tau.1.h.f.aux(t)                    + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t)                   + omega*Tau.3.c.f.aux(t))*r1.normal.h.c.exp                        - lambda*c( Tau.2.delayed.h.f(t), Tau.2.delayed.c.f(t) )*r1.normal.h.c.exp                                  - c(alpha.h, alpha.c)*r1.normal.h.c.exp                                                                                                                                                               + c(epsilon.1.h, epsilon.1.c)*S.normal.h.c.exp)
            dr1.amplificado.h.c      =   c(1,1)*(1/r1.amplificado.h.c.exp)*(c(theta.h.f(t)*r1.amplificado.h.c.exp[2]             - mu.discharge*r1.amplificado.h.c.exp[1],        - theta.h.f(t)*r1.amplificado.h.c.exp[2]       + theta.c(t) + mu.discharge*r1.amplificado.h.c.exp[1] )         - mu.death*r1.amplificado.h.c.exp                                                                                                                                                                                                                                                                                                                          + c(Tau.1.h.f.aux(t)                    + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t)                   + omega*Tau.3.c.f.aux(t))*kleb.prob.func(t)*r1.normal.h.c.exp                                                                                                                  - c(alpha.h, alpha.c)*r1.amplificado.h.c.exp                                                                                                                                                                                                                                                         - lambda*(c(1,1)/(1/365.25+tau.2) - c(phi.h, phi.c))*r1.amplificado.h.c.exp  - c(phi.h, phi.c)*r1.amplificado.h.c.exp)
            dr1.amplificado.beta.h.c =   c(1,1)*(1/r1.amplificado.beta.h.c.exp)*(c(theta.h.f(t)*r1.amplificado.beta.h.c.exp[2]   - mu.discharge*r1.amplificado.beta.h.c.exp[1],   - theta.h.f(t)*r1.amplificado.beta.h.c.exp[2]  + theta.c(t) + mu.discharge*r1.amplificado.beta.h.c.exp[1] )    - mu.death*r1.amplificado.beta.h.c.exp  + (1/c(N[1],N[2]))* c(beta.1.h, beta.1.c)*(r1.normal.h.c.exp + nu.inf*r1.amplificado.h.c.exp + nu.inf*r1.amplificado.beta.h.c.exp)*nu.susc*S.amplificado.h.c.exp                                                                                                                                                   + c(Tau.1.h.f.aux(t)                    + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t)                   + omega*Tau.3.c.f.aux(t))*(1-kleb.prob.func(t))*r1.normal.h.c.exp  - lambda*c( Tau.2.delayed.h.f(t), Tau.2.delayed.c.f(t) )*r1.amplificado.beta.h.c.exp                        - c(alpha.h, alpha.c)*r1.amplificado.beta.h.c.exp                                                                                                                                                     + c(epsilon.1.h, epsilon.1.c)*nu.susc*S.amplificado.h.c.exp                                                                                                                 + c(phi.h, phi.c)*r1.amplificado.h.c.exp)

            ## 
            dr2.normal.h.c           =   c(1,1)*(1/r2.normal.h.c.exp)*(c(theta.h.f(t)*r2.normal.h.c.exp[2]                       - mu.discharge*r2.normal.h.c.exp[1],             - theta.h.f(t)*r2.normal.h.c.exp[2]            + theta.c(t) + mu.discharge*r2.normal.h.c.exp[1] )              - mu.death*r2.normal.h.c.exp            + (1/c(N[1],N[2]))* c(beta.2.h, beta.2.c)*(r2.normal.h.c.exp + nu.inf*(r2.amplificado.h.c.exp + r2.amplificado.beta.h.c.exp) )*S.normal.h.c.exp                                                                                                                                                                    - c(Tau.1.h.f.aux(t) + Tau.2.f.r(t)[1]  + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t) + Tau.2.f.r(t)[2] + omega*Tau.3.c.f.aux(t))*r2.normal.h.c.exp                                                                                                                                    - c(alpha.h, alpha.c)*r2.normal.h.c.exp                                                                                                                                                               + c(epsilon.2.h, epsilon.2.c)*S.normal.h.c.exp)
            dr2.amplificado.h.c      =   c(1,1)*(1/r2.amplificado.h.c.exp)*(c(theta.h.f(t)*r2.amplificado.h.c.exp[2]             - mu.discharge*r2.amplificado.h.c.exp[1],        - theta.h.f(t)*r2.amplificado.h.c.exp[2]       + theta.c(t) + mu.discharge*r2.amplificado.h.c.exp[1] )         - mu.death*r2.amplificado.h.c.exp                                                                                                                                                                                                                                                                                                                          + c(Tau.1.h.f.aux(t) + Tau.2.f.r(t)[1]  + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t) + Tau.2.f.r(t)[2] + omega*Tau.3.c.f.aux(t))*kleb.prob.func(t)*r2.normal.h.c.exp                                                                                                                  - c(alpha.h, alpha.c)*r2.amplificado.h.c.exp)
            dr2.amplificado.beta.h.c =   c(1,1)*(1/r2.amplificado.beta.h.c.exp)*(c(theta.h.f(t)*r2.amplificado.beta.h.c.exp[2]   - mu.discharge*r2.amplificado.beta.h.c.exp[1],   - theta.h.f(t)*r2.amplificado.beta.h.c.exp[2]  + theta.c(t) + mu.discharge*r2.amplificado.beta.h.c.exp[1] )    - mu.death*r2.amplificado.beta.h.c.exp                                                                                                                                                     + (1/c(N[1],N[2]))*c(beta.2.h, beta.2.c)*(r2.normal.h.c.exp + nu.inf*(r2.amplificado.h.c.exp + r2.amplificado.beta.h.c.exp))*nu.susc*S.amplificado.h.c.exp      + c(Tau.1.h.f.aux(t) + Tau.2.f.r(t)[1]  + omega*Tau.3.h.f.aux(t), Tau.1.c.f.aux(t) + Tau.2.f.r(t)[2] + omega*Tau.3.c.f.aux(t))*(1-kleb.prob.func(t))*r2.normal.h.c.exp                                                                                                              - c(alpha.h, alpha.c)*r2.amplificado.beta.h.c.exp                                                                                                                                                     + c(epsilon.2.h, epsilon.2.c)*nu.susc*S.amplificado.h.c.exp)                                 

            ## estimate solutions
            sol = list(c(dS.normal.h.c, dS.amplificado.h.c, dr1.normal.h.c, dr1.amplificado.h.c, dr1.amplificado.beta.h.c, dr2.normal.h.c, dr2.amplificado.h.c, dr2.amplificado.beta.h.c),
                       N = N, N.total                    = sum(N),
                       prop.S.normal.h.c.exp             = prop.S.normal.h.c.exp,
                       prop.S.amplificado.h.c.exp        = prop.S.amplificado.h.c.exp,
                       prop.r1.normal.h.c.exp            = prop.r1.normal.h.c.exp,
                       prop.r1.amplificado.h.c.exp       = prop.r1.amplificado.h.c.exp,
                       prop.r1.amplificado.beta.h.c.exp  = prop.r1.amplificado.beta.h.c.exp,
                       prop.r2.normal.h.c.exp            = prop.r2.normal.h.c.exp,
                       prop.r2.amplificado.h.c.exp       = prop.r2.amplificado.h.c.exp,
                       prop.r2.amplificado.beta.h.c.exp  = prop.r2.amplificado.beta.h.c.exp,
                       prop.r1.h                         = prop.r1.h.c[1],
                       prop.r1.c                         = prop.r1.h.c[2],
                       kleb.restricted.percent           = kleb.restricted.percent
                       )
            
        })}
        
        
        state =  c(S.normal.h.c  = with(parameters0,S.normal.init),  S.amplificado.h.c = with(parameters0,S.amplificado.init),
                   r1.normal.h.c = with(parameters0,r1.normal.init), r1.amplificado.h.c = with(parameters0,r1.amplificado.init), r1.amplificado.beta.h.c = with(parameters0,r1.amplificado.beta.init),
                   r2.normal.h.c = with(parameters0,r2.normal.init), r2.amplificado.h.c = with(parameters0,r2.amplificado.init), r2.amplificado.beta.h.c = with(parameters0,r2.amplificado.beta.init)
                   )
        
        ll = 1/12                                ## time step --> monthly
        times = seq(initime, endtime, by = ll)
        out =  ode(y = state, times = times, func = diffEq, parms = parameters)
        as.data.frame(out)
    }

    ## define costs functions
    AMR.klebsiella.cost =  function(parameters.candidates.ob){
        parameters0[names(parameters.candidates.ob)]   =  parameters.candidates.ob
        
        paux = as.numeric(parameters0[names(parameters.candidates.ob)])
        names(paux) = names(parameters0[names(parameters.candidates.ob)])
        
        if (length(paux) == 3){ paux = c(paux, epsilon.slope.to.equivalent.f(paux["epsilon.1.c.slope"])); names(paux)[3] = "epsilon.1.c.equivalent"}
        message("parms in costs: ");print(t(as.matrix(paux)))
        
        out     = bacterial.master(parameters =  parameters0)
        costs   = modCost(model = out, obs =  prevalence.resistance.ceftriaxon.dir.outpatient[,c("time", "prop.r1.c", "std")], err =  "std")
        costs.1 = modCost(model = out, obs =  prevalence.resistance.ceftriaxon.dir.inpatient[,c("time", "prop.r1.h", "std")],  err =  "std", cost = costs)
        costs.1
    }
    
    AMR.klebsiella.cost.full =  function(parameters.candidates.ob){
        parameters0[names(parameters.candidates.ob)]   =  parameters.candidates.ob
        out     = bacterial.master(parameters =  parameters0)
        costs   = modCost(model = out, obs =  prevalence.resistance.ceftriaxon.dir.outpatient[,c("time", "prop.r1.c", "std")], err =  "std")
        costs.1 = modCost(model = out, obs =  prevalence.resistance.ceftriaxon.dir.inpatient[,c("time", "prop.r1.h", "std")],  err =  "std", cost = costs)
        costs.1
    }
    
    AMR.klebsiella.cost.2 = function(lpars) {out = AMR.klebsiella.cost(exp(lpars)); message("\n COSTS$model: ", out$model, "\n");out}

    if(quiere.sensitivity == TRUE){ 
        message("For generalized sensitivity")
    }
    
    
    if(init.solution == TRUE){
        message("init Solutions will be computed")
        sol      = bacterial.master(parameters = parameters0)
        costs.0 =  AMR.klebsiella.cost.full(parameters.candidates.ob = parameters0[c("beta.1.h", "beta.1.c")])
        message("costs.0$model"); print(costs.0$model)
        
        pp.r1.h = ggplot(prevalence.resistance.ceftriaxon.dir.inpatient, aes(x = time, y = prop.r1.h)) + geom_point(col = "red", size = 2) + ggtitle("Hospital") +
            geom_line(data = sol, aes(x =  time, y = prop.r1.h), size = 1.25, col =  "orange", lty = 2) + ylab("Fraction") + xlab("") +
            geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.inpatient, aes(x =  time,  ymin = ci.lower, ymax= ci.upper), colour =  "red")
        
        pp.r1.c = ggplot(prevalence.resistance.ceftriaxon.dir.outpatient, aes(x = time, y = prop.r1.c)) + geom_point(col = "red", size = 2) + ggtitle("Community") +
            geom_line(data = sol, aes(x =  time, y = prop.r1.c), size = 1.25, col =  "orange", lty = 2) +
            geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.outpatient, aes(x =  time,  ymin = ci.lower, ymax= ci.upper), colour =  "red") +
            ylab("Fraction") + xlab("Years since 2000")

        x11(title = "init solution");grid.arrange(pp.r1.h,pp.r1.c, nrow = 1)
        fn.aux = paste(fitfolder,"init_sol.tiff", sep = "/")
        tiff(fn.aux, height = 5, width = 5, unit = "in", res = 300);print(grid.arrange(pp.r1.h,pp.r1.c, nrow = 1));dev.off()        
        write.table(sol, "init_solution")
    }

    ## optimizsation boundaries
    lower. = c("beta.1.h" = (-10),   "beta.1.c" = (-10),   "epsilon.1.c.slope" = (-20))
    upper. = c("beta.1.h" = log(30), "beta.1.c" = log(30), "epsilon.1.c.slope" = log(-2*log(1-r1.last.data*1)/(t.r1.last.data*t.r1.last.data)))

    ## set to given value if to run a single epsilon, or as dummy to name files from optimizations 
    if((optimize.epsilon) + (baseline.epsilon) == 1)  epsilons.share = epsilon.equivalent.exogeno   

    ## main fits
    if(quiere.fits == TRUE){
        message("Yes, quiere.fits")
        if(optimize.epsilon == TRUE){
            fit.pars = matrix(NA, nrow = length(epsilons.share) , ncol =  3);
            colnames(fit.pars) = c("ln(beta.1.h)",   "ln(beta.1.c)", "ln(epsilon.1.c.slope)")
        }else{
            fit.pars = matrix(NA, nrow = length(epsilons.share) , ncol =  2);
            colnames(fit.pars) = c("ln(beta.1.h)",   "ln(beta.1.c)")
        }
        
        sols.list = list()
        
        for(i in 1:length(epsilons.share)){
            message("iteration, i= ", i  , " of ", length(epsilons.share))
            parameters0["epsilon.1.c.slope"] = -2*log(1-r1.last.data*epsilons.share[i])/(t.r1.last.data*t.r1.last.data)
            message("Dios mio, va a optimizar!")
            endtime = 17.5
            
            if (optimize.epsilon == TRUE) parameters.candidates["epsilon.1.c.slope"] = parameters0["epsilon.1.c.slope"]
            
            print("parameters.candidates before fitting")
            print(parameters.candidates)
            if (optimize.epsilon == TRUE) message("epsilon equivalent: ", epsilon.slope.to.equivalent.f(as.numeric(parameters.candidates["epsilon.1.c.slope"])))
            print("end fitting parameters")

            Sys.time.init = Sys.time()
            
            Fit = modFit(p =  log(unlist(parameters.candidates)),
                         f =  AMR.klebsiella.cost.2,
                         method = "L-BFGS-B", lower = lower.[names(parameters.candidates)],  upper = upper.[names(parameters.candidates)]
                         )
            
            Fit.all.outcomes  = Fit
            
            print("summary(Fit)");print(summary(Fit))
            
            Sys.time.final =  Sys.time()
            lasted = Sys.time.final - Sys.time.init
            message("Fit lasted: "); print(lasted)
            
            estimate = exp(summary(Fit)$par[,"Estimate"])
            ## the parameters with ci:
            ci.parameters.95 = data.frame(estimate, lower =  estimate*exp(-1.96*summary(Fit)$par[,"Std. Error"]), upper =  estimate*exp(1.96*summary(Fit)$par[,"Std. Error"]))
            print(" ci.parameters.95");print( ci.parameters.95)
            
            write.table(Fit$par,                 paste(fitfolder,paste("Fit_par",                 100*epsilons.share[i], sep =  "_"), sep = "/"))
            write.table(Fit$var_ms_unweighted,   paste(fitfolder,paste("Fit_var_ms_unweighted",   100*epsilons.share[i], sep =  "_"), sep = "/"))
            write.table(ci.parameters.95,        paste(fitfolder,paste("Fit_ci_par",              100*epsilons.share[i], sep =  "_"), sep = "/"))
            write.table(summary(Fit)$cov.scaled, paste(fitfolder,paste("summary(Fit)_cov_scaled", 100*epsilons.share[i], sep =  "_"), sep = "/"))
            write.table(summary(Fit)$par,        paste(fitfolder,paste("summary(Fit)_par_std",    100*epsilons.share[i], sep =  "_"), sep = "/"))
            message("Optimization finalized")
        }
    }## closes fits


    ## MCMC parameters sampling 
    if(mcmc.compute == TRUE){
        
        for(i in 1:length(epsilons.share)){
            
            message("iteration, i= ", i  , " of ", length(epsilons.share))
            
            parameters0["epsilon.1.c.slope"] = -2*log(1-r1.last.data*epsilons.share[i])/(t.r1.last.data*t.r1.last.data)
            
            var0. <-          read.table(paste(fitfolder, paste("Fit_var_ms_unweighted",   100*epsilons.share[i], sep =  "_"), sep = "/"))
            cov0 <- as.matrix(read.table(paste(fitfolder, paste("summary(Fit)_cov_scaled", 100*epsilons.share[i], sep =  "_"), sep = "/"))* 2.4^2/3)
            Fit.  <-          read.table(paste(fitfolder, paste("Fit_par",                 100*epsilons.share[i], sep =  "_"), sep = "/"))
            Fit =  Fit.$x
            names(Fit) = rownames(Fit.)
            var0  = var0.$x
            names(var0) = rownames(var0.)      

            Sys.time.init = Sys.time()
            
            try(        MCMC <- modMCMC(f = AMR.klebsiella.cost.2,
                                        p = Fit,
                                        niter = 1000,
                                        jump = cov0,
                                        var0 = var0,
                                        wvar0 = 1,
                                        updatecov = 25,
                                        lower = (lower.)[names(parameters.candidates)],
                                        upper = (upper.)[names(parameters.candidates)]
                                        ), silent = TRUE
                )
            
            Sys.time.final = Sys.time()
            message("MCMC generation took:");  print(Sys.time.final - Sys.time.init)
            
            
            if( length(grep("MCMC", ls()))  == 0 ) {
                warning("MCMC failed, MCMC$pars will be set to optimization point estimates")
                MCMC =  list( pars = data.frame(t(Fit.)),  bestpar =  as.data.frame((Fit.)))      
                print("MCMC");print(MCMC);print(class(MCMC))
                write.table(Sys.time(), "MCMC_fail.log")
            }else{
                MCMC.init =  MCMC
                print(summary(MCMC))
                plot(MCMC, Full = TRUE)
                pdf(paste(fitfolder, "MCMC_chain.pdf",sep = "/")); plot(MCMC, Full = TRUE);dev.off()
            }    
            MCMC$pars = exp(MCMC$pars)                                
            write.table(MCMC$pars,    paste(fitfolder, paste("MCMC_pars", 100*epsilons.share[i], sep =  "_"), sep  = "/"))
            write.table(exp(MCMC$bestpar), paste(fitfolder, paste("MCMC_bestpars", 100*epsilons.share[i], sep =  "_"), sep  = "/"))
        }
    }## closes MCMC sampling
    
    ## solve differential equations with the mcmc parameters
    if(mcmc.solve == TRUE){
        
        if(amc.scenarios == FALSE){factors.amc = 1}else{     factors.amc           =  c(1+0.5, 1+0.25, 1+0.1, 1, 1-0.1, 1-0.25, 1-0.5)}
        if(beta.1.h.scenarios == FALSE){beta.1.h.factors = 1}else{beta.1.h.factors =  c(1+0.5, 1+0.25, 1+0.1, 1, 1-0.1, 1-0.25, 1-0.5)}
        
        mcmcRange = function(func = bacterial.master, parms = parameters.candidates, parInput){
            parameters = parameters0
            ii =  1:dim(parInput)[1]
            by2dataframe.f(lapply(ii, function(j) {
                parameters[names(parameters.candidates)] = parInput[j,names(parameters.candidates)];
                print(parameters[names(parameters.candidates)])
                out = func(parameters)[,c("time", "prop.r1.h", "prop.r1.c")];
                message("iteration: ", j, " of: ", length(ii))
                out$iteration = j
                out$type      = parInput[j,"type"]
                out
            }))}
        
        factors.amc = expand.grid(factors.amc, beta.1.h.factors); colnames(factors.amc) = c("factors.amc","beta.1.h.factor")
        print("factors.amc");print(factors.amc)
        
        for(i in 1:length(epsilons.share)){
            parameters0["epsilon.1.c.slope"] = epsilon.equivalent.to.slope.f(epsilons.share[i])
            MCMC.pars =  unique(read.table(paste(fitfolder, paste("MCMC_pars", 100*epsilons.share[i], sep =  "_"), sep  = "/")))
            print("Accepted chains: "); print(dim(MCMC.pars))
            fit.pars       = exp(read.table(paste(fitfolder, paste("Fit_par",                 100*epsilons.share[i], sep =  "_"), sep = "/")))
            MCMC.pars$type = 0
            bestline       = as.data.frame(cbind(t(fit.pars), type = 1) )
            MCMC.pars      = rbind(MCMC.pars, bestline)
            ## write.table(MCMC.pars, "MCMC_pars_withbest")
            
            endtime = 25
            
            message("solutions iteration begins")
            solutions.mcmc =     apply(factors.amc, 1, function(fi){
                print("amc.factor:fi treatment scenario"); print(fi)
                
                bacterial.master.aux =  function(parameters = parameters.iteration, amc.factor = fi["factors.amc"], beta.1.h.factor = fi["beta.1.h.factor"], amc.sc.onlyrestricted =  amc.sc.onlyrestricted.st) bacterial.master(parameters = parameters, amc.factor = amc.factor, beta.1.h.factor = beta.1.h.factor, amc.sc.onlyrestricted = amc.sc.onlyrestricted )
                solution = mcmcRange(func = bacterial.master.aux, parms = parameters.candidates, parInput = MCMC.pars)
                
                solution$factor.amc   = fi[1]
                solution$beta.1.h.amc = fi[2]
                solution
            }
            )
            
            solutions.mcmc.df = by2dataframe.f(solutions.mcmc)
            solutions.mcmc.df.melt = melt(solutions.mcmc.df, id.vars =  c("time",  "iteration", "type", "factor.amc", "beta.1.h.amc"))
            
            solutions.mcmc.df.melt$place[grep("\\.c\\.|\\.c$", solutions.mcmc.df.melt$variable)] = "Community"
            solutions.mcmc.df.melt$place[grep("\\.h\\.|\\.h$", solutions.mcmc.df.melt$variable)] = "Hospitals"
            
            write.table(solutions.mcmc.df.melt, paste(fitfolder, paste(figure.label, 100*epsilons.share[i], sep =  "_"), sep  = "/"))
            
            rm(solutions.mcmc.df.melt)
        }
    } ## closes mcmc solutions


    ## plots(paper figures)
    if (mcmc.plot == TRUE){
        prevalence.resistance.ceftriaxon.dir.inpatient =  read.table("prevalence_resistance_ceftriaxon_dir_inpatient")
        prevalence.resistance.ceftriaxon.dir.outpatient = read.table("prevalence_resistance_ceftriaxon_dir_outpatient")
        prevalence.resistance.ceftriaxon.dir.inpatient$place =  "Hospitals"
        prevalence.resistance.ceftriaxon.dir.outpatient$place   =  "Community"
        
        if(amc.scenarios == FALSE){factors.amc = 1}else{     factors.amc           =  c(1+0.5, 1+0.25, 1+0.1, 1, 1-0.1, 1-0.25, 1-0.5)}
        if(beta.1.h.scenarios == FALSE){beta.1.h.factors = 1}else{beta.1.h.factors =  c(1+0.5, 1+0.25, 1+0.1, 1, 1-0.1, 1-0.25, 1-0.5)}
        
        factors.amc = expand.grid(factors.amc, beta.1.h.factors); colnames(factors.amc) = c("factors.amc","beta.1.h.factor")
        
        if(figure.label ==  "Figure2"){
            
            for(i in 1:length(epsilons.share)){
                aux.sim =                     read.table(paste(fitfolder, paste(figure.label, 100*epsilons.share[i], sep =  "_"), sep  = "/"))
                
            aux.sim$epsilon.equivalent =  100*epsilons.share[i]

            if (i == 1){range.sim = aux.sim}else{range.sim = rbind(range.sim, aux.sim)}
            }
            
            data.main.ci = range.sim;
            rm(range.sim)

            ## Adding the simulations with optimized epsilon
            
            data.main.ci.optim = subset(read.table( paste("FiguresG3_nuinf", "_", nu.inf.sens.parallel.value, "/","Figure3A_50",sep = "")  ), factor.amc ==1)
            aux.fit.optim = exp(t(read.table( paste("FiguresG3_nuinf", "_", nu.inf.sens.parallel.value, "/","Fit_par_50", sep = "")  ) )[1, "epsilon.1.c.slope"])

            data.main.ci.optim$epsilon.equivalent = 100*epsilon.slope.to.equivalent.f(aux.fit.optim)
            
            data.main.ci =  rbind(data.main.ci, data.main.ci.optim)
            data.main.ci =  data.main.ci[order(data.main.ci[,"epsilon.equivalent"]),]
            
            data.main.ci = subset(data.main.ci, iteration == 1) ##plot only the optimal value
            
        }else{ ## figure.label !=  "Figure2"        
            data.main.ci = read.table(paste(fitfolder, paste(figure.label, 100*epsilons.share, sep =  "_"), sep  = "/"))
            data.main.ci$epsilon.equivalent =  100*epsilons.share 
        }
        
        epsilon.equivalent =  100*epsilons.share
        
        aux.f =  as.factor(data.main.ci$beta.1.h.amc)
        data.main.ci$beta.1.h.amc = aux.f
        aux.f =  as.factor(data.main.ci$factor.amc)
        data.main.ci$factor.amc = aux.f
        aux.f =  as.factor(data.main.ci$place)
        data.main.ci$place = aux.f
        aux.f =  as.factor(data.main.ci$iteration)
        data.main.ci$iteration = aux.f
        
        data.main.ci$beta.1.h.amc = factor(data.main.ci$beta.1.h.amc, levels =  c("1.5", "1.25", "1.1","1", "0.9", "0.75", "0.5"))
        data.main.ci$factor.amc   = factor(data.main.ci$factor.amc,   levels =  c("1.5", "1.25", "1.1","1", "0.9", "0.75", "0.5"))
        
        pp.lines = ggplot(data.main.ci) +
            facet_grid(as.factor(place)~., scales = "free") +
            stat_summary(data = data.main.ci, aes(x =  as.numeric(time) + 2000, y = as.numeric(100*value), colour =  as.factor(epsilon.equivalent)), geom="point",   fun.y    =  median, size = 1  ) +
            geom_point(data = prevalence.resistance.ceftriaxon.dir.inpatient,     aes(x = time + 2000,   y =  100*prop.r1.h),                     colour =  "grey") +
            geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.inpatient,  aes(x = time + 2000,  ymin = 100*ci.lower, ymax= 100*ci.upper), colour =  "grey") +
            geom_point(data = prevalence.resistance.ceftriaxon.dir.outpatient,    aes(x = time + 2000,   y =  100*prop.r1.c),                     colour =  "grey") +
            geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.outpatient, aes(x = time + 2000,  ymin = 100*ci.lower, ymax= 100*ci.upper), colour =  "grey") +                              
            theme(
                legend.position = "right",
                legend.key.width = unit(2, "line"),
                legend.text =  element_text(size = 15),
                legend.title =  element_text(size = 15),
                plot.title = element_text(lineheight=.8, face="bold", size = 18),
                axis.line.x = element_line(colour = "black"),
                axis.line.y = element_line(colour = "black"),
                axis.line = element_line(colour = "black"),
                panel.background = element_blank(),
                axis.text.x = element_text(colour =  "black", size = 13, hjust = 1, angle = 90),
                axis.text.y = element_text(colour =  "black", size = 17.5),
                axis.title.x = element_text(colour =  "black", size = 20),
                axis.title.y = element_text(colour =  "black", size = 20),
                strip.text.x = element_text(size = 20, colour = "black")
            ) + coord_cartesian(ylim =  ylimite, xlim = c(2000,2025)) +
            ylab("Prevalence of Resistance (%)") + xlab("Year")
        
        x11(title = "epsilon sensitivity") ; print(pp.lines)
        
        
        if (figure.label != "Figure2"){
            
            if(amc.scenarios == TRUE)     {     labels.amc = c("Increased by 50%", "Increased by 25%", "Increased by 10%", "Stable", "Reduced by 10%", "Reduced by 25%", "Reduced by 50%")}else {labels.amc = "Stable"}
            if(beta.1.h.scenarios == TRUE){labels.beta.1.h = c("Increased by 50%", "Increased by 25%", "Increased by 10%", "Stable", "Reduced by 10%", "Reduced by 25%", "Reduced by 50%")}else{labels.beta.1.h = "Stable"}
            
            if(beta.1.h.scenarios == TRUE){
                data.aux.beta = data.main.ci
                pp.lines.2.0 = ggplot(data = data.aux.beta ) +
                    geom_point(aes(x =  time + 2000, y = 100*value, pch = as.factor(beta.1.h.amc) ), col  = topo.colors(7)[4], size = 0.5) + scale_linetype_manual(name = "Scenario of \n in hospital transmission by 2025",  labels = labels.beta.1.h, values =  c(2,3,4,1,4,3,2) )               
            }else{
                pp.lines.2.0 = ggplot(data = data.main.ci, aes(colour =  as.factor(factor.amc))) +
                    geom_point(aes(x =  time + 2000, y = 100*value, lty = beta.1.h.amc), size = 0.5) +
                    scale_colour_manual(name = "Scenario of \n antimicrobial consumption by 2025", labels = labels.amc, values =  (topo.colors(length(levels(as.factor(factors.amc$factors.amc))))))
            }
            
            pp.lines.2 = pp.lines.2.0 +
                facet_grid(~as.factor(place), scales = "free") +
                geom_point(data = prevalence.resistance.ceftriaxon.dir.inpatient, aes(x = time + 2000,   y =  100*prop.r1.h),  colour =  "grey", size = 0.1, alpha = 0.5)   +
                geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.inpatient, aes(x =  time + 2000,  ymin = 100*ci.lower, ymax= 100*ci.upper), colour =  "grey")     +
                geom_point(data = prevalence.resistance.ceftriaxon.dir.outpatient, aes(x = time + 2000,   y =  100*prop.r1.c),                           colour =  "grey")  +
                geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.outpatient, aes(x =  time + 2000,  ymin = 100*ci.lower, ymax= 100*ci.upper), colour =  "grey")    +
                theme(
                    legend.position = "right",
                    legend.key.width = unit(2, "line"),
                    legend.text =  element_text(size = 15),
                    legend.title =  element_text(size = 15),
                    plot.title = element_text(lineheight=.8, face="bold", size = 18),
                    axis.line.x = element_line(colour = "black"),
                    axis.line.y = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    panel.background = element_blank(),
                    axis.text.x = element_text(colour =  "black", size = 13, hjust = 1, angle = 90),
                    axis.text.y = element_text(colour =  "black", size = 17.5),
                    axis.title.x = element_text(colour =  "black", size = 20),
                    axis.title.y = element_text(colour =  "black", size = 20),
                    strip.text.x = element_text(size = 20, colour = "black")
                ) + ylab("Prevalence of Resistance (%)") + xlab("Year") +
                coord_cartesian(ylim =  ylimite, xlim = NULL)
            
            alpha = ifelse(amc.sc.onlyrestricted.st ==1, 0.1, 0.3)
            
            if(beta.1.h.scenarios == FALSE){  ## scenarios other than hospital transmission
                pp.ci =
                    ggplot(data =  data.main.ci,  size = 0.5) +
                    facet_grid(~as.factor(place), scales = "free") +
                    geom_line(data = subset(data.main.ci, type ==1), aes(x =  time + 2000, y = 100*value, lty = beta.1.h.amc, colour =  as.factor(factor.amc)), size = 1  ) +
                    scale_fill_manual(name = "Scenario of \n antimicrobial consumption by 2025", labels = labels.amc, values =  (topo.colors(length(levels(as.factor(factors.amc$factors.amc)))))) +
                    scale_colour_manual(name = "Scenario of \n antimicrobial consumption by 2025", labels = labels.amc, values =  (topo.colors(length(levels(as.factor(factors.amc$factors.amc)))))) +
                    geom_point(data = prevalence.resistance.ceftriaxon.dir.inpatient,     aes(x = time + 2000,     y = 100*prop.r1.h), colour =  "darkgrey")   +
                    geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.inpatient,  aes(x = time + 2000,  ymin = 100*ci.lower,   ymax= 100*ci.upper), colour =  "darkgrey") +
                    geom_point(data = prevalence.resistance.ceftriaxon.dir.outpatient,    aes(x = time + 2000,     y = 100*prop.r1.c), colour =  "darkgrey")  +
                    geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.outpatient, aes(x = time + 2000,  ymin = 100*ci.lower,   ymax= 100*ci.upper), colour =  "darkgrey") +
                    theme(
                        legend.position = "bootom",
                        legend.key.width = unit(2, "line"),
                        legend.text =  element_text(size = 15),
                        legend.title =  element_text(size = 15),
                        plot.title = element_text(lineheight=.8, face="bold", size = 18),
                        axis.line.x = element_line(colour = "black"),
                        axis.line.y = element_line(colour = "black"),
                        axis.line = element_line(colour = "black"),
                        panel.background = element_blank(),
                        axis.text.x = element_text(colour =  "black", size = 13, hjust = 1, angle = 90),
                        axis.text.y = element_text(colour =  "black", size = 17.5),
                        axis.title.x = element_text(colour =  "black", size = 20),
                        axis.title.y = element_text(colour =  "black", size = 20),
                        strip.text.x = element_text(size = 20, colour = "black")
                    ) + ylab("Prevalence of Resistance (%)") + xlab("Year") +
                    coord_cartesian(ylim =  ylimite, xlim = NULL) + ylab("Prevalence of Resistance (%)") + xlab("Year")
                pp.ci
                pp.ci.zoom = pp.ci + coord_cartesian(ylim =  c(2, 16), xlim = c(2019,2025))
                pp.ci.zoom                
            }else{    ## hospital transmission scenarios
                alpha = 0.15
                
                pp.ci =
                    ggplot(data =  data.main.ci,  size = 0.5) +
                    facet_grid(~as.factor(place), scales = "free") +                    
                    geom_line(data = subset(data.main.ci, type ==1), aes(x =  time + 2000, y = 100*value, lty = beta.1.h.amc), colour =  topo.colors(7)[4], size = 1  ) +
                    scale_fill_manual(name = "Scenario of \n antimicrobial consumption by 2025", labels = labels.amc, values =  (topo.colors(length(levels(as.factor(factors.amc$factors.amc)))))) +
                    scale_colour_manual(name = "Scenario of \n antimicrobial consumption by 2025", labels = labels.amc, values =  (topo.colors(length(levels(as.factor(factors.amc$factors.amc)))))) +
                    geom_point(data = prevalence.resistance.ceftriaxon.dir.inpatient,     aes(x = time + 2000,     y = 100*prop.r1.h), colour =  "darkgrey")   +
                    geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.inpatient,  aes(x = time + 2000,  ymin = 100*ci.lower,   ymax= 100*ci.upper), colour =  "darkgrey") +
                    geom_point(data = prevalence.resistance.ceftriaxon.dir.outpatient,    aes(x = time + 2000,     y = 100*prop.r1.c), colour =  "darkgrey")  +
                    geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.outpatient, aes(x = time + 2000,  ymin = 100*ci.lower,   ymax= 100*ci.upper), colour =  "darkgrey") +
                    theme(
                        legend.position = "none",
                        legend.key.width = unit(2, "line"),
                        legend.text =  element_text(size = 15),
                        legend.title =  element_text(size = 15),
                        plot.title = element_text(lineheight=.8, face="bold", size = 18),
                        axis.line.x = element_line(colour = "black"),
                        axis.line.y = element_line(colour = "black"),
                        axis.line = element_line(colour = "black"),
                        panel.background = element_blank(),
                        axis.text.x = element_text(colour =  "black", size = 13, hjust = 1, angle = 90),
                        axis.text.y = element_text(colour =  "black", size = 17.5),
                        axis.title.x = element_text(colour =  "black", size = 20),
                        axis.title.y = element_text(colour =  "black", size = 20),
                        strip.text.x = element_text(size = 20, colour = "black")
                    ) + ylab("Prevalence of Resistance (%)") + xlab("Year") +
                    coord_cartesian(ylim =  ylimite, xlim = NULL) + ylab("Prevalence of Resistance (%)") + xlab("Year")  +
                    scale_linetype_manual(name = "Scenario of \n in hospital transmission by 2025",  labels = labels.beta.1.h, values =  c(2,3,4,1,4,3,2) )
                pp.ci
                pp.ci.zoom = pp.ci + coord_cartesian(ylim =  c(2, 16), xlim = c(2019,2025))
                pp.ci.zoom
            }
            x11(title = "lines prevalence mcmc"); print(pp.ci)
            
            tiff(tiffname1, height = 5, width = 5, unit = "in", res = 300);print(pp.ci);dev.off()
            tiff(tiffname2, height = 5, width = 5, unit = "in", res = 300);print(pp.ci.zoom);dev.off()
            
            data.aux.2019 =  subset(data.main.ci, time == 19)
            data.aux.2025 =  subset(data.main.ci, time == 25)
            
            data.aux.2019.2025 =  data.aux.2019
            data.aux.2019.2025$increase =  (data.aux.2025$value - data.aux.2019$value)/data.aux.2019$value
            data.aux.2019.2025$factor.amc.sc[data.aux.2019.2025$factor.amc == 1.5] = "Increased by 50%"
            data.aux.2019.2025$factor.amc.sc[data.aux.2019.2025$factor.amc == 1.25] = "Increased by 25%"
            data.aux.2019.2025$factor.amc.sc[data.aux.2019.2025$factor.amc == 1.1] = "Increased by 10%"
            data.aux.2019.2025$factor.amc.sc[data.aux.2019.2025$factor.amc == 1] = "Stable"
            data.aux.2019.2025$factor.amc.sc[data.aux.2019.2025$factor.amc == 0.9] = "Reduced by 10%"
            data.aux.2019.2025$factor.amc.sc[data.aux.2019.2025$factor.amc == 0.75] = "Reduced by 25%"
            data.aux.2019.2025$factor.amc.sc[data.aux.2019.2025$factor.amc == 0.5] = "Reduced by 50%"
            
            ## --summary outcome from the terminal data:
            data.aux.2019.2025.point.estimates.from.best.probability = subset(data.aux.2019.2025, (type==1))
            print("data.aux.2019.2025.point.estimates.from.best.probability"); print(data.aux.2019.2025.point.estimates.from.best.probability)
            ## confidence intervales from ther resampled mcmc
            increases.ci95 =  by2dataframe.f(by(data.aux.2019.2025, data.aux.2019.2025[,c("factor.amc", "place")], function(x) data.frame(factor.amc = x[1,"factor.amc"], place = x[1,"place"], quantile(x$increase, c(0.025,0.5,0.975)))))
            colnames(increases.ci95)[3] = "increases.ci95"
            print("increases.ci95"); print(increases.ci95)

            if(beta.1.h.scenarios == TRUE){
                increases.ci95 =  by2dataframe.f(by(data.aux.2019.2025, data.aux.2019.2025[,c("beta.1.h.amc", "place")], function(x) data.frame(factor.amc = x[1,"beta.1.h.amc"], place = x[1,"place"], quantile(x$increase, c(0.025,0.5,0.975)))))
                colnames(increases.ci95)[3] = "increases.ci95"
                print(increases.ci95)
            }

            # calculate error bars from mcmc for plotting
            data.aux.subset = subset(data.aux.2019.2025)
            intervals = by2dataframe.f(by(data.aux.subset,
                                          data.aux.subset[,c("factor.amc", "place")],
                                          function(x) data.frame(factor.amc = x[1,"factor.amc"],
                                                                 place = x[1, "place"],
                                                                 min = quantile(x$increase,0.025),
                                                                 median = quantile(x$increase,0.50),
                                                                 max = quantile(x$increase,0.975)) 
                                          ))
            
            ## bars with 2025 comparison values
            
            if(beta.1.h.scenarios == FALSE){  ## scenarios other than hospital transmission

                pp.terminal = ggplot(intervals, aes(x = as.factor(factor.amc), y = median,  fill = as.factor(factor.amc) ) ) +
                    facet_grid(~as.factor(place), scales = "free") +
                    geom_bar(aes(x=as.factor(factor.amc), y=median), lty = 1, position=position_dodge(width=0.95), stat =  "identity") +
                    geom_errorbar(aes(ymin=min, ymax=max), width=.2, lwd = 0.5, position=position_dodge(.9)) +
                    scale_fill_manual(name = "Scenario of \n antimicrobial consumption by 2025", labels = labels.amc, values =  (topo.colors(length(levels(as.factor(factors.amc$factors.amc)))))) +
                    theme(
                        legend.position = "none",
                        legend.key.width = unit(2, "line"),
                        legend.text =  element_text(size = 15),
                        legend.title =  element_text(size = 15),
                        plot.title = element_text(lineheight=.8, face="bold", size = 18),
                        axis.line.x = element_line(colour = "black"),
                        axis.line.y = element_line(colour = "black"),
                        axis.line = element_line(colour = "black"),
                        panel.background = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(colour =  "black", size = 10),
                        axis.title.x = element_text(colour =  "black", size = 20),
                        axis.title.y = element_text(colour =  "black", size = 20),
                        strip.text.x = element_text(size = 20, colour = "black")
                    ) + ylab("Increase 2019-2025(%)") + xlab("")
                
                if(figure.label == "Figure3A") {pp.terminal =  pp.terminal + scale_y_continuous(limits=c(-1,13.3), breaks=  c(seq(-1, 3,by =0.5), seq(4,13.3,by = 2)))}
                pp.terminal
            }else{ ## hospital transmission scenarios
                pp.terminal = ggplot(data.aux.2019.2025, aes(x = (as.factor(beta.1.h.amc)) , y = as.numeric(increase),  fill = as.factor(factor.amc) ) ) +
                    facet_grid(~as.factor(place), scales = "free") +
                    geom_bar(data = subset(data.aux.2019.2025, type ==1), lty = 1, position=position_dodge(width=0.95), stat =  "identity", fill =  "black") +
                    stat_summary(fun.y =  median, geom="bar", position=position_dodge(width=0.95), fill =  "grey", alpha = 0.5) +
                    stat_summary(fun.data =  function(x) median_hilow(x, conf.int=.95), geom = "errorbar", position=position_dodge(width=0.95), size = 0.75, col = topo.colors(7)[4], aes(lty = as.factor(beta.1.h.amc))) +
                    scale_fill_manual(name = "Scenario of \n antimicrobial consumption by 2025", labels = labels.amc, values =  (topo.colors(length(levels(as.factor(factors.amc$factors.amc)))))) +
                    theme(
                        legend.position = "right",
                        legend.key.width = unit(2, "line"),
                        legend.text =  element_text(size = 15),
                        legend.title =  element_text(size = 15),
                        plot.title = element_text(lineheight=.8, face="bold", size = 18),
                        axis.line.x = element_line(colour = "black"),
                        axis.line.y = element_line(colour = "black"),
                        axis.line = element_line(colour = "black"),
                        panel.background = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(colour =  "black", size = 13),
                        axis.title.x = element_text(colour =  "black", size = 20),
                        axis.title.y = element_text(colour =  "black", size = 20),
                        strip.text.x = element_text(size = 20, colour = "black")
                    ) + ylab("Increase 2019-2025(%)") + xlab("") + scale_y_continuous(limits=c(-0.61, 1.7)) + scale_linetype_manual(name = "Scenario of \n in hospital transmission by 2025",  labels = labels.beta.1.h, values =  c(2,3,4,1,4,3,2) )
                
                pp.terminal
            }
            
            x11(title = "pp terminal"); print(pp.terminal)
            tiff(paste(fitfolder,paste(figure.label, "terminal.tiff", sep = "_"), sep = "/"), height = 5, width = 5, unit = "in", res = 300);print(pp.terminal);dev.off()
            
        }else{ ## figure.label == "Figure2"
            
            data.main.ci. = data.main.ci
            data.main.ci.$optim = 0
            data.main.ci.$optim[ data.main.ci.$epsilon.equivalent == setdiff(data.main.ci.$epsilon.equivalent,100*epsilons.share)] = 1
            
            pp.ci =
                ggplot(data =  data.main.ci.,  size = 0.5) +
                facet_grid(~as.factor(place), scales = "free") +
                stat_summary(data = subset(data.main.ci., time>=0), aes(x =  time + 2000, y = 100*value, colour =  as.factor(epsilon.equivalent), lty = as.factor(optim)),                                geom="line",   fun.y    =  median, size = 1  ) +
                geom_point(data = prevalence.resistance.ceftriaxon.dir.inpatient,     aes(x = time + 2000,     y = 100*prop.r1.h), colour =  "darkgrey")   +
                geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.inpatient,  aes(x = time + 2000,  ymin = 100*ci.lower,   ymax= 100*ci.upper), colour =  "darkgrey") +
                geom_point(data = prevalence.resistance.ceftriaxon.dir.outpatient,    aes(x = time + 2000,     y = 100*prop.r1.c), colour =  "darkgrey")  +
                geom_errorbar(data = prevalence.resistance.ceftriaxon.dir.outpatient, aes(x = time + 2000,  ymin = 100*ci.lower,   ymax= 100*ci.upper), colour =  "darkgrey") +
                theme(
                    legend.position = "none",
                    legend.key.width = unit(2, "line"),
                    legend.text =  element_text(size = 15),
                    legend.title =  element_text(size = 15),
                    plot.title = element_text(lineheight=.8, face="bold", size = 18),
                    axis.line.x = element_line(colour = "black"),
                    axis.line.y = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    panel.background = element_blank(),
                    axis.text.x = element_text(colour =  "black", size = 13, hjust = 1, angle = 90),
                    axis.text.y = element_text(colour =  "black", size = 17.5),
                    axis.title.x = element_text(colour =  "black", size = 20),
                    axis.title.y = element_text(colour =  "black", size = 20),
                    strip.text.x = element_text(size = 20, colour = "black")
                ) + ylab("Prevalence of Resistance (%)") + xlab("Year") +
                coord_cartesian(ylim =  ylimite, xlim = NULL) + ylab("Prevalence of Resistance (%)") + xlab("Year") +
                scale_colour_manual(name = "External force \n of colonization equivalent*", labels =  round(unique(data.main.ci.$epsilon.equivalent),2),
                                    values =  rev(rainbow(length(levels(as.factor(data.main.ci.$epsilon.equivalent)))))) +
                scale_linetype_manual(name = " ", labels =  c("Assumption", "Model Fit"), values =  c(1,2))            
            x11(title = "pp ci"); print(pp.ci)
            
            tiff(tiffname1, height = 5, width = 5, unit = "in", res = 300);print(pp.ci);dev.off()

            data.aux.2019 =  subset(data.main.ci, time == 19)
            data.aux.2025 =  subset(data.main.ci, time == 25)
            
            data.aux.2019.2025 =  data.aux.2019
            data.aux.2019.2025$increase =  (data.aux.2025$value - data.aux.2019$value)/data.aux.2019$value
            
            pp.terminal = ggplot(data.aux.2019.2025, aes(x = as.factor(epsilon.equivalent) , y = as.numeric(increase),  fill = as.factor(epsilon.equivalent) ) ) +
                facet_grid(~as.factor(place), scales = "free") +
                stat_summary(fun.y =  median, geom="bar", lty = 1, position=position_dodge(width=0.95)) +
                stat_summary(fun.data =  function(x) median_hilow(x, conf.int=.95), geom = "errorbar", lty = 1, position=position_dodge(width=0.95), size = 0.25, col = "black") +
                theme(
                    legend.position = "right",
                    legend.key.width = unit(2, "line"),
                    legend.text =  element_text(size = 15),
                    legend.title =  element_text(size = 15),
                    plot.title = element_text(lineheight=.8, face="bold", size = 18),
                    axis.line.x = element_line(colour = "black"),
                    axis.line.y = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    panel.background = element_blank(),
                    axis.text.x = element_blank(),
                    axis.text.y = element_text(colour =  "black", size = 13),
                    axis.title.x = element_text(colour =  "black", size = 20),
                    axis.title.y = element_text(colour =  "black", size = 20),
                    strip.text.x = element_text(size = 20, colour = "black")
                ) + ylab("Increase 2019-2025(%)") + xlab("") 
            
            x11(title = "terminal"); print(pp.terminal)
            
            tiff(paste(fitfolder,paste(figure.label, "terminal.tiff", sep = "_"), sep = "/"), height = 5, width = 5, unit = "in", res = 300);print(pp.terminal);dev.off()
        }
        

        data.aux.extrems =  subset(data.main.ci, ( (time ==5.5) + (time ==19) + (time ==25) + (time ==17.5)) ==1 )
        
        extrems.ci = by2dataframe.f(by(data.aux.extrems, data.aux.extrems[,c("time", "factor.amc", "beta.1.h.amc",  "variable",  "place", "epsilon.equivalent")], function(x) {
            aux1 = x[1,c("time", "factor.amc", "beta.1.h.amc",  "variable",  "place", "epsilon.equivalent")]
            aux2 = t(as.data.frame(quantile(x$value, c(0,0.025,0.25,0.5,0.75,0.975,1))))
            cbind(aux1,aux2)
        } ))
        
        write.table(extrems.ci, paste(fitfolder, paste(figure.label, "extrem_ci", sep = "_"), sep =  "/") )
        
    }#closes mcmc.plot
    
}#closes the function for parallel execution
