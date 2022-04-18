##--Figures s8,s9 
rm(list = ls())

require("parallel")
require(snowfall)

start.time = Sys.time()

run.the.fit = TRUE
run.mcmc = TRUE
solve.mcmc.scenarios = TRUE

if(run.the.fit == TRUE){
    sfInit(parallel = TRUE, cpus = 6)
    sfLapply(c(0, seq(0.45,0.65, by =0.05)),  function(esp, nuifnsens = 2){
        fitfolder <<- paste("FiguresG3_nuinf_", nuifnsens, "/eps_",100*esp, sep = "")
        source("core.r");
        parallel.fits.f(nu.inf.sensitivity.logical = FALSE,
                        nu.inf.sens.parallel.value = nuifnsens,
                        epsilon.equivalent.sensitivity = esp,
                        init.solution = FALSE,
                        quiere.sensitivity = FALSE,
                        fit.incl.epsilon = FALSE,
                        fit.excl.epsilon = TRUE,
                        baseline.epsilon = TRUE,
                        mcmc.compute =FALSE,
                        mcmc.solve = FALSE,
                        mcmc.plot = FALSE,
                        optimize.epsilon = FALSE,
                        fitfolder = fitfolder
                        )
    } )
    sfStop()
}

if(run.mcmc == TRUE){
    sfInit(parallel = TRUE, cpus = 6)    
    sfLapply(c(0, seq(0.45,0.65, by =0.05)),  function(esp, nuifnsens = 2){
        fitfolder <<- paste("FiguresG3_nuinf_", nuifnsens, "/eps_",100*esp, sep = "")
        source("core.r");        
        parallel.fits.f(nu.inf.sensitivity.logical = FALSE,
                        nu.inf.sens.parallel.value = nuifnsens,
                        epsilon.equivalent.sensitivity = esp,
                        init.solution = FALSE,
                        quiere.sensitivity = FALSE,
                        fit.incl.epsilon = FALSE,
                        fit.excl.epsilon = FALSE,
                        baseline.epsilon = TRUE,
                        mcmc.compute =TRUE,
                        mcmc.solve = FALSE,
                        mcmc.plot = FALSE,
                        optimize.epsilon = FALSE,
                        fitfolder = fitfolder
                        )
    } )
    sfStop()
}

if(solve.mcmc.scenarios == TRUE){
    sfInit(parallel = TRUE, cpus = 6)    
    sfLapply(c(0, seq(0.45,0.65, by =0.05)),  function(esp, nuifnsens = 2){
        fitfolder <<- paste("FiguresG3_nuinf_", nuifnsens, "/eps_",100*esp, sep = "")       
        source("core.r")                
        parallel.fits.f(nu.inf.sensitivity.logical = FALSE,
                        nu.inf.sens.parallel.value = nuifnsens,
                        epsilon.equivalent.sensitivity = esp,
                        init.solution = FALSE,
                        quiere.sensitivity = FALSE,
                        fit.incl.epsilon = FALSE,
                        fit.excl.epsilon = FALSE,
                        baseline.epsilon = TRUE,
                        mcmc.compute = FALSE,
                        mcmc.solve = TRUE,
                        mcmc.plot = TRUE,
                        optimize.epsilon = FALSE,
                        fitfolder = fitfolder,
                        Figure = "3A"
                        )
    } )
    sfStop()
}

end.time = Sys.time()
execution.time = end.time - start.time

message("execution.time: ", execution.time)
