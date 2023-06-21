using VAM
using DataFrames
using Distributions
m=5
dfcov=DataFrame(cov1=rand(Uniform(),m), cov2=rand(Uniform(),m), cov3=rand(Uniform(),m))

m = @vam( time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5| 1*cov1 + -2cov2 + 3cov3)))
df = rand(m, 29, datacov=dfcov)

ml = mle(m,df,dfcov)
contrast(ml.mle)
params(ml.mle.model)
m.params_cov
m.vars_cov
m.datacov