
rm(list = ls())

require("parallel")
require(snowfall)

start.time = Sys.time()


nuifnsens = 2
fitfolder <<- paste("testcode", nuifnsens, sep = "_")
source("core.r") 
parallel.fits.f(nu.inf.sens.parallel.value = nuifnsens,
                epsilon.equivalent.sensitivity = 0.5,
                nu.inf.sensitivity.logical = TRUE,
                init.solution = TRUE,
                quiere.sensitivity =FALSE,
                fit.incl.epsilon = FALSE,
                baseline.epsilon = FALSE,
                mcmc.compute = FALSE,
                mcmc.solve = FALSE,
                mcmc.plot = FALSE,
                optimize.epsilon = TRUE,
                fitfolder = fitfolder
                )

end.time = Sys.time()
execution.time = end.time - start.time
message("execution.time: ", execution.time)


aa = read.table("init_solution")

pp = ggplot(data = aa) + geom_line(aes(x = time, y = ))
