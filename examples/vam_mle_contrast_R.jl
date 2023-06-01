using VAM
using DataFrames
using RCall
using Test
using OrderedCollections

mutable struct ModelResult 
	models::Dict
	results::Dict
	ModelResult() = new(Dict(),Dict())
end

modres = ModelResult()

modres.models = Dict(
	:WInf => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))),
		:r => [
			"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
			"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
			] 
	),
	:W1 => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))),
		:r => [
			"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
			"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
			] 
	)
)


function update!(modres::ModelResult, key::Symbol)
	model = modres.models[key]
	data = model[:data]
	θ=model[:θ]
	result = Dict()
	m = VAM.MLE(model[:vam], data)
	R"""
	require(VAM)
	simData <- eval(parse(text=$(model[:r][1])))
	form <-  eval(parse(text=$(model[:r][2])))
	mle <- mle.vam(form,data=simData)
	theta <- $(model[:θ])
	res <- list()
	res$lnL <- logLik(mle,theta,TRUE,FALSE,FALSE)
	res$dlnL<- logLik(mle,theta,FALSE,TRUE,FALSE)
	res$d2lnL <- logLik(mle,theta,FALSE,FALSE,TRUE)
	res$C <- contrast(mle,theta,TRUE,FALSE,FALSE)
	res$dC <- contrast(mle,theta,FALSE,TRUE,FALSE)
	res$d2C <- contrast(mle,theta,FALSE,FALSE,TRUE)
	"""
	result[:resR] = @rget res
	result[:resJL] = OrderedDict(
		:lnL => contrast(m, θ, alpha_fixed=true),
		:dlnL => gradient(m, θ, alpha_fixed=true),
		:d2lnL => hessian(m, θ, alpha_fixed=true),
		:C => contrast(m, θ),
		:dC => gradient(m, θ)[2:end],
		:d2C => hessian(m, θ)[2:end,2:end]
	)
	modres.results[key]=result
end

function test(modres::ModelResult)
	for (key, result) in modres.results
		println("model: $key")
		for k in [:lnL, :dlnL, :d2lnL, :C, :dC]
			@test result[:resJL][k] == result[:resR][k]
		end
	end
end

update!(modres, :W1)
test(modres)