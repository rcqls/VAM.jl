using VAM
using DataFrames
using Optim
using Random
using FreqTables
using CSV

Random.seed!(3)

## RCpp: sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
m = @vam(Temps & Type ~ (ARA1(.9) | Weibull(0.01,2.5)) & (ARAInf(.7) | AtIntensity(0.2)))
df = simulate(m, @stop(s < 1000))
freqtable(df.type)
params(m)
#@run mle(m, params(m), df; method=Newton())
res = mle(m, params(m), df; method=Newton())
Optim.minimizer(res)
res = mle(m, params(m), df, method=LBFGS())
Optim.minimizer(res)

contrast(m, params(m), df)
gradient(m, params(m), df)
hessian(m, params(m), df)