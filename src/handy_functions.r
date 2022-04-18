##handy functions along the algorithm. List/by type to data.frame and integral/derivative
by2dataframe.f =  function(by.ob) {for(i in 1:length(by.ob)){if(i == 1) {by.ob.df = by.ob[[1]]}else{by.ob.df = rbind(by.ob.df, by.ob[[i]])}}; by.ob.df  }
epsilon.slope.to.equivalent.f = function(ep.slope){ round((1-exp(-(t.r1.last.data*t.r1.last.data/2)*ep.slope))/r1.last.data,3) }
epsilon.equivalent.to.slope.f =  function(ep.equivalent) -2*log(1-r1.last.data*ep.equivalent)/(t.r1.last.data*t.r1.last.data)
