require(VAM)
simCMPM<-sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.7) | AtIntensity(0.2)))
simData<-simulate(simCMPM, Size>=1000)
print(head(simData))
print(table(simData$Type))
mleCMPM <- mle.vam( Time & Type ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.7) | AtIntensity(0.2)), data=simData)
print(coef(mleCMPM))
print(mleCMPM$fixed)
theta<-c(0.001,2.5,0.9, .7)
print(logLik(mleCMPM,theta,FALSE,TRUE,FALSE))
print(contrast(mleCMPM,theta,FALSE,TRUE,FALSE))