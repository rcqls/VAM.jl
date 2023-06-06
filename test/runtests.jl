using Test
using VAM
using DataFrames
using RCall
using OrderedCollections

include("testing.jl")

modtest = ModelTest()

insert!(modtest, 
	:W_ARA∞ => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	),
	:W_ARA∞_ARA∞ => Dict(
		:θ =>  [0.3,1.8,0.3,0.8],
		:data => DataFrame(Temps=[3.36, 4.1, 6.1, 7.2],Type=[-1, 1, 1, -1]),
		:vam => @vam(Temps & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)) & (ARA∞(0.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	), 
	:W_ARA1 => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
		# 	"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
		# ] 
	),
	:W_ARA∞bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)))
	),
	:W_ARA1bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)))
	),
	:LL_ARA∞ => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA∞(0.4) | LogLinear(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	), 
	:LL_ARA1 => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | LogLinear(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
		# 	"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
		# ] 
	),
	:LL_ARA∞bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | LogLinear(0.001,2.5)))
	),
	:LL_ARA1bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA1(0.4) | LogLinear(0.001,2.5)))
	),
	:W_ARA∞_ARA∞_ARA1 => Dict(
		:θ => [0.3,1.8,0.3,0.4,0.7],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16,7.16],Type=[1,2,2,1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)) & (ARA∞(0.5)+ARA1(0.7)))
	),
	:W_ARA∞_AGAN_ARA1 => Dict(
		:θ => [0.3,1.8,0.3,0.7],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16,7.16],Type=[1,2,2,1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)) & (AGAN()+ARA1(0.7)))
	),
    :W_ARA∞_MS => Dict(
        :θ => [0.3,0.8,0.6],
        :data => DataFrame(System=[1,2],Time=[3.36,2.34],Type=[-1,-1]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))) 
    ),
    :W_ARA∞_ARA∞_MS => Dict(
        :θ => [0.3,1.8,0.3,0.8],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5))) 
    )
)

modtest.models[:W_ARA∞_ARA∞_ARA1][:r][2]

 

#empty!(modtest)
update!(modtest, :W_ARA∞_ARA∞_ARA1)
update!(modtest, :W_ARA∞_ARA∞)
update!(modtest)
test(modtest)

modtest.results[:W_ARA∞_ARA∞_MS][:r]
modtest.results[:W_ARA∞_ARA∞_MS][:jl]
# string(modtest.models[:W1][:vam].formula)
# rterms(modtest.models[:W1][:vam], DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]))