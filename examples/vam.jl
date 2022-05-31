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
θ = [0.3,1.8,0.6]
lnL = -2.30449245951301
dlnL = [-5.52619699756427,-1.45367181592636,0]
d2lnL = vcat([[-11.1111111111111,-10.7372278181901,0],[-10.7372278181901,-4.21250787723964,0],[0,0,0]]...)
c = -1.62415430907299
dC = [0.555555555555556,0]
d2C = reshape([-0.308641975308642,0,0,0],2,2)
VAM.contrast(m,θ)
VAM.gradient(m, θ)
VAM.hessian(m, θ)
VAM.contrast(m, θ, alpha_fixed=true)
VAM.gradient(m, θ, alpha_fixed=true)
VAM.hessian(m, θ, alpha_fixed=true)
@run VAM.contrast(m,θ)
@run VAM.gradient(m, θ)
@run VAM.hessian(m, θ)

data = DataFrame(time=[3.36],type=[-1])
m = mle(
	@vam(Time & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5))), 
	data
)
θ = [0.3,0.8,0.6]
lnL = -3.65431355894635
dlnL = [-13.7944691820681, -8.74189899224908, 0]
d2lnL = vcat([[-11.1111111111111,-40.3396633074969,0],[-40.3396633074969,-31.9886643027399, 0], [0,0,0]])
C = -1.15270302129339
dC = [1.00478465516967,0]
d2C = reshape([-0.678445876029819,0,0,0], 2, 2)
VAM.contrast(m,θ)
VAM.gradient(m, θ)
VAM.hessian(m, θ)
VAM.contrast(m, θ, alpha_fixed=true)
VAM.gradient(m, θ, alpha_fixed=true)
VAM.hessian(m, θ, alpha_fixed=true)

data = DataFrame(time=[3.36],type=[-1])
m = mle(
	@vam(Time & Type ~ (ARAInf(0.4) | Weibull3(0.001,2.5,5.0))),
	data
)
θ = [0.3,1.8,4,0.6]
lnL =  -6.28349650594271
dlnL = [-20.8805277090384,-14.1662361334174,-0.920550636729322,0]
d2lnL = vcat([[-11.1111111111111,-55.726172072379,-3.43082096301078,0],[-55.726172072379,-36.7534932861936,-3.48854152770058,0],[-3.43082096301078,-3.48854152770058,0.0228198059801746],[0,0,0,0,0]])
C = -2.00229062844351
dC = [0.250199288093325,-0.0329926542472398,0]
dC = vcat([[-0.0292037862530474,-0.0369910685383343,0],[-0.0369910685383343,0.01048162449677,0],[0,0,0]])
VAM.contrast(m,θ)
VAM.gradient(m, θ)
VAM.hessian(m, θ)
VAM.contrast(m, θ, alpha_fixed=true)
VAM.gradient(m, θ, alpha_fixed=true)
VAM.hessian(m, θ, alpha_fixed=true)


data = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = mle(
	@vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))), 
	data
)
θ = [0.3,0.8,0.6]
lnL = -7.37830963135462
dlnL = [9.33348796771948,5.77076155284033,1.19923836457015]
d2lnL = vcat([[-44.4444444444444,-5.26230480903023,-0.723044448779292],[-5.26230480903023,-7.7471885008781,-6.30726171420585],[-0.723044448779292,-6.30726171420585,0.435684125802398]])
C = -5.36231016699152
dC = [2.08694474533596,0.693079297460398]
d2C = vcat([-4.31732300351524,-3.55104259186756],[-3.55104259186756,-0.351517905180765])
VAM.contrast(m,θ)
VAM.gradient(m, θ)
VAM.hessian(m, θ)
VAM.contrast(m, θ, alpha_fixed=true)
VAM.gradient(m, θ, alpha_fixed=true)
VAM.hessian(m, θ, alpha_fixed=true)


data = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = mle(
	@vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))), 
	data
)
θ = [0.3,0.8,0.6]
lnL = -7.60531410020218
dlnL = [9.43046924122556,6.95299388770745,0.950641033424141]
d2lnL = vcat([[-44.4444444444444,-5.58984454112956,-0.51705781427391],[-5.58984454112956,-8.18607969766148,-5.04449078731948],[-0.51705781427391,-5.04449078731948,1.70955451742894]])
C = -5.5202288755881
dC = [2.90098061507367,0.575831839583211]
d2C = vcat([[-4.65895384634265,-3.11529381117442],[-3.11529381117442,0.845549840945982]])
VAM.contrast(m,θ)
VAM.gradient(m, θ)
VAM.hessian(m, θ)
VAM.contrast(m, θ, alpha_fixed=true)
VAM.gradient(m, θ, alpha_fixed=true)
VAM.hessian(m, θ, alpha_fixed=true)


data = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = mle(
	@vam(Time & Type ~ (AGAN() | Weibull(0.001,2.5))), 
	data
)
θ = [0.3,0.8]
lnL = -6.90098251895707
dlnL = [8.75359407434948,3.3717767832448]
d2lnL = vcat([[-44.4444444444444,-2.40399936753948],[-2.40399936753948,-7.6652713149061]])
C = -5.25256034408912
dC = 1.99329429144225
d2C = -9.26821752088532
VAM.contrast(m,θ)
VAM.gradient(m, θ)
VAM.hessian(m, θ)
VAM.contrast(m, θ, alpha_fixed=true)
VAM.gradient(m, θ, alpha_fixed=true)
VAM.hessian(m, θ, alpha_fixed=true)


data = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = mle(
	@vam(Time & Type ~ (ABAO() | LogLinear(0.001,2.5))), 
	data
)
θ = [0.3,0.8]
lnL = -13.6870254780068
dlnL = [-62.9837808690102,-73.9249749593491]
d2lnL = vcat([[-44.4444444444444,-304.849916531164],[-304.849916531164,-390.943849373404]])
C = -1.77041141396439
dC = 1.55193690275558
d2C = -4.47702275378921
VAM.contrast(m,θ)
VAM.gradient(m, θ)
VAM.hessian(m, θ)
VAM.contrast(m, θ, alpha_fixed=true)
VAM.gradient(m, θ, alpha_fixed=true)
VAM.hessian(m, θ, alpha_fixed=true)
