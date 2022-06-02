using VAM
using DataFrames
using Optim
m = @vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)))
df = simulate(m, @stop(s < 100))
params(m)
res = mle(m, params(m), df)
Optim.minimizer(res)
