using VAM

m = @vam( time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5| 1*cov1 + -2cov2 + 3cov3)))
m.params_cov
m.vars_cov