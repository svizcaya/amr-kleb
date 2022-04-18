##--Figures 3, s5, s6
rm(list = ls())

require("parallel")
require(snowfall)

start.time = Sys.time()

run.the.fit = TRUE
run.mcmc = TRUE
solve.mcmc.scenarios = TRUE

if(run.the.fit == TRUE){
    sfInit(parallel = TRUE, cpus = 3)
    sfLapply(seq(1, 3, by = 1), function(nuifnsens){
        fitfolder <<- paste("FiguresG3_nuinf", nuifnsens, sep = "_")
        source("core.r") 
        parallel.fits.f(nu.inf.sens.parallel.value = nuifnsens,
                        epsilon.equivalent.sensitivity = 0.5,
                        nu.inf.sensitivity.logical = TRUE,
                        init.solution = FALSE,
                        quiere.sensitivity =FALSE,
                        fit.incl.epsilon = TRUE,
                        baseline.epsilon = FALSE,
                        mcmc.compute = FALSE,
                        mcmc.solve = FALSE,
                        mcmc.plot = FALSE,
                        optimize.epsilon = TRUE,
                        fitfolder = fitfolder
                        )
    } )
    sfStop()
}


if(run.mcmc == TRUE){
    sfInit(parallel = TRUE, cpus = 3)
    sfLapply(seq(1, 3, by = 1), function(nuifnsens){
        
        fitfolder <<- paste("FiguresG3_nuinf", nuifnsens, sep = "_")
        source("core.r") 
        parallel.fits.f(nu.inf.sens.parallel.value = nuifnsens,
                        epsilon.equivalent.sensitivity = 0.5,
                        nu.inf.sensitivity.logical = TRUE,
                        init.solution = FALSE,
                        quiere.sensitivity =FALSE,
                        fit.incl.epsilon = FALSE,
                        baseline.epsilon = FALSE,
                        mcmc.compute = TRUE,
                        mcmc.solve = FALSE,
                        mcmc.plot = FALSE,
                        optimize.epsilon = TRUE,
                        fitfolder = fitfolder
                        )
    } )
    sfStop()
}

if(solve.mcmc.scenarios == TRUE){
    sfInit(parallel = TRUE, cpus = 3)
    
    sfLapply(seq(1, 3, by = 1), function(nuifnsens){        
        fitfolder <<- paste("FiguresG3_nuinf", nuifnsens, sep = "_")
        source("core.r") 
        parallel.fits.f(
            nu.inf.sens.parallel.value = nuifnsens,
            epsilon.equivalent.sensitivity = 0.5,
            nu.inf.sensitivity.logical = TRUE,
            init.solution = FALSE,
            quiere.sensitivity =FALSE,
            fit.incl.epsilon = FALSE,
            baseline.epsilon = FALSE,
            mcmc.compute = FALSE,
            mcmc.solve = TRUE,
            mcmc.plot = TRUE,
            optimize.epsilon = TRUE,
            fitfolder = fitfolder,
            Figure = "3C")
    } )
    sfStop()
}

end.time = Sys.time()
execution.time = end.time - start.time
message("execution.time: ", execution.time)
