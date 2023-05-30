using VAM
using DataFrames
using RCall

data = DataFrame(time=[3.36],type=[-1])
m = VAM.MLE(
	@vam(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5))), 
	data
)
θ = [0.3,1.8,0.6]
lnL = -2.30449245951301
dlnL = [-5.52619699756427,-1.45367181592636,0]
d2lnL = [
	-11.1111111111111 -10.7372278181901 0;
	-10.7372278181901 -4.21250787723964 0;
	0 0 0
]
lnLJ = contrast(m, θ, alpha_fixed=true)
gradient(m, θ, alpha_fixed=true)
hessian(m, θ, alpha_fixed=true)
c = -1.62415430907299
dC = [0.555555555555556,0]
d2C = [-0.308641975308642 0; 0 0]
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)

R"""
require(VAM)
simData<-data.frame(Time=c(3.36),Type=c(-1),row.names=1:1)
mle <- mle.vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=simData)
theta<-c(0.3,1.8,0.6)
lnL<- -2.30449245951301
dlnL<- c(-5.52619699756427,-1.45367181592636,0)
d2lnL<- matrix(c(-11.1111111111111,-10.7372278181901,0,-10.7372278181901,-4.21250787723964,0,0,0,0),nrow=3,byrow=TRUE)
C<- -1.62415430907299
dC<-c(0.555555555555556,0)
d2C<-matrix(c(-0.308641975308642,0,0,0),nrow=2,byrow=TRUE)
lnLR <- logLik(mle,theta,TRUE,FALSE,FALSE) #,equals(lnL,tolerance=0.00000000000001))
dlnLR<- logLik(mle,theta,FALSE,TRUE,FALSE) #,equals(dlnL,tolerance=0.00000000000001))
d2lnLR <- logLik(mle,theta,FALSE,FALSE,TRUE)#,equals(d2lnL,tolerance=0.00000000000001))
CR <- contrast(mle,theta,TRUE,FALSE,FALSE)#,equals(C,tolerance=0.00000000000001))
dCR <- contrast(mle,theta,FALSE,TRUE,FALSE)#,equals(dC,tolerance=0.00000000000001))
d2CR <- contrast(mle,theta,FALSE,FALSE,TRUE)#,equals(d2C,tolerance=0.00000000000001))
"""

@rget simData
@rget lnLR