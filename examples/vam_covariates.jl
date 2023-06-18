using VAM
using DataFrames
using Distributions
m=5
df_cov=DataFrame(cov1=rand(Uniform(),m), cov2=rand(Uniform(),m), cov3=rand(Uniform(),m))

m = @vam( time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5| 1*cov1 + -2cov2 + 3cov3)))
rand(m, 29, data_cov=df_cov)


VAM.covariates!(m, df_cov)
m.params_cov
m.vars_cov
m.data_cov