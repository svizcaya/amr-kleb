## template example for steering core routine
rm(list = ls())
source("core.r");
parallel.fits.f(nu.inf.sensitivity.logical = TRUE,
                nu.inf.sens.parallel.value = 2,
                epsilon.equivalent.sensitivity = 0,
                init.solution = FALSE,
                quiere.sensitivity = FALSE,
                fit.incl.epsilon = FALSE,
                fit.excl.epsilon = FALSE,
                baseline.epsilon = TRUE,
                mcmc.compute = FALSE,
                mcmc.solve = TRUE,
                mcmc.plot = FALSE,  
                optimize.epsilon = FALSE,
                fitfolder = "Results",
                Figure = "2"
                )
} )
