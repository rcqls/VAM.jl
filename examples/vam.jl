using VAM
using DataFrames

s = simulator(
	@vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2))),
  	@stop (size < 30) & (time < 30)
)

df = simulate(s)

params(s.model)

df2 = simulate(s,system=30)

m = mle(
	@vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2))),
	df2
)
m.model.k
m.model.data

m.model.time

select_data(m.model,2)

m.model.time
data = DataFrame(time=[3.36],type=[-1])
m = mle(
	@vam(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5))), 
	data
)
m.model.k
θ = [0.3,1.8,0.6]
lnL = -2.30449245951301
dlnL = [-5.52619699756427,-1.45367181592636,0]
d2lnL = vcat([[-11.1111111111111,-10.7372278181901,0],[-10.7372278181901,-4.21250787723964,0],[0,0,0]]...)
c = -1.62415430907299
dC = [0.555555555555556,0]
d2C = reshape([-0.308641975308642,0,0,0],2,2)
# @test logLik(mle,theta,TRUE,FALSE,FALSE),equals(lnL,tolerance=0.00000000000001))
# expect_that(logLik(mle,theta,FALSE,TRUE,FALSE),equals(dlnL,tolerance=0.00000000000001))
# expect_that(logLik(mle,theta,FALSE,FALSE,TRUE),equals(d2lnL,tolerance=0.00000000000001))
VAM.contrast(m,θ)
VAM.contrast(m,θ, alpha_fixed=true)
VAM.gradient(m, θ)
@run VAM.contrast(m,θ)
@run VAM.gradient(m, θ)
params(m.model)