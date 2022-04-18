##--Figure s4A
rm(list = ls())

require("parallel")
require(snowfall)
start.time = Sys.time()

run.the.fit = TRUE
run.mcmc = TRUE
solve.mcmc.scenarios = TRUE
plot.scenarios.epsilon = TRUE

if(run.the.fit == TRUE){
    sfInit(parallel = TRUE, cpus = 5)
    sfLapply(seq(0, 0.6, by = 0.15),  function(esp, nuifnsens = 2){
        fitfolder <<- paste("FiguresG2_nuinf", nuifnsens, sep = "_")
        source("core.r");
        parallel.fits.f(nu.inf.sensitivity.logical = TRUE,
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
    sfInit(parallel = TRUE, cpus = 5)
    sfLapply(seq(0, 0.6, by = 0.15),  function(esp, nuifnsens = 2){
        fitfolder <<- paste("FiguresG2_nuinf", nuifnsens, sep = "_")
        source("core.r");
        parallel.fits.f(nu.inf.sensitivity.logical = TRUE,
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

if(solve.mcmc.scenarios == TRUE) {
    sfInit(parallel = TRUE, cpus = 5)
    sfLapply(seq(0, 0.6, by = 0.15),  function(esp, nuifnsens = 2){
        fitfolder <<- paste("FiguresG2_nuinf", nuifnsens, sep = "_")
        source("core.r");
        parallel.fits.f(nu.inf.sensitivity.logical = TRUE,
                        nu.inf.sens.parallel.value = nuifnsens,
                        epsilon.equivalent.sensitivity = esp,
                        init.solution = FALSE,
                        quiere.sensitivity = FALSE,
                        fit.incl.epsilon = FALSE,
                        fit.excl.epsilon = FALSE,
                        baseline.epsilon = TRUE,
                        mcmc.compute = FALSE,
                        mcmc.solve = TRUE,
                        mcmc.plot = FALSE,  
                        optimize.epsilon = FALSE,
                        fitfolder = fitfolder,
                        Figure = "2"
                        )
    } )
    sfStop()    
}

if( plot.scenarios.epsilon == TRUE) {
    nuifnsens = 2
    fitfolder <<- paste("FiguresG2_nuinf", nuifnsens, sep = "_")
    source("core.r");
    parallel.fits.f(nu.inf.sensitivity.logical = TRUE,
                    nu.inf.sens.parallel.value = nuifnsens,
                    epsilon.equivalent.sensitivity = 0.5, 
                    init.solution = FALSE,
                    quiere.sensitivity = FALSE,
                    fit.incl.epsilon = FALSE,
                    fit.excl.epsilon = FALSE,
                    baseline.epsilon = FALSE,
                    mcmc.compute = FALSE,
                    mcmc.solve = FALSE,
                    mcmc.plot = TRUE, 
                    optimize.epsilon = FALSE,
                    fitfolder = fitfolder,
                    Figure = "2"
                    )   
}

end.time = Sys.time()
execution.time = end.time - start.time
message("execution.time: ", execution.time)
