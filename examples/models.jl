using VAM

vam = @vam(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5)))
VAM.init!(vam)

ex_f= :(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) )
ex_f2 = :(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2)))
ex_f.head
ex_f.args
ex_f.args[3].head
ex_f.args[3].args
ex_f2.args[3].args
VAM.parse_model(ex_f)
VAM.parse_model(ex_f2)