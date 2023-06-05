using VAM
using DataFrames
using RData
using RCall



R"""
require(VAM)
data("AMC_Amb")
"""

@rget AMC_Amb

R"""
AMC_mle<-mle.vam(Time & Type ~ (ARAInf(0.6)|Weibull(1,3)),data=AMC_Amb)
res <- coef(AMC_mle)
cAMC <- contrast(AMC_mle,res,TRUE,FALSE,FALSE)
gAMC <- contrast(AMC_mle,res,FALSE,TRUE,FALSE)
"""
@rget res
@rget cAMC
@rget gAMC
m = @vam(Time & Type ~ (ARAInf(0.6) | Weibull(1.0,3.0)))
mle(m, AMC_Amb)
params(m)
