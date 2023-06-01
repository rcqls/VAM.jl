using VAM
using DataFrames
using RCall
using Test
using OrderedCollections

const ModelDict = Dict{Symbol, Any}
function empty!(md::ModelDict)
	for (k, _) in md
		delete!(md, k)
	end
end



function rterms(m::VAM.Model, data::DataFrame)
	f = m.formula
	df_rexpr = string("data.frame(", 
		join(map(names(data)) do var
			string(var,"=c(",join(data[!,var],","),")")
		end,","),
	")")
	vam_repxr = string(replace(string(f.args[2]),"(" => "", ")" => ""), " ~ (", replace(string(f.args[3]),"FamilyModel" => ""), ")")
	[df_rexpr, vam_repxr]
end
mutable struct ModelTest
	models::ModelDict
	results::ModelDict
	ModelTest() = new(ModelDict(),ModelDict())
end

modtest = ModelTest()

function insert!(modtest::ModelTest, models::Vararg{Pair{Symbol, ModelDict}})
	for model in models
		m = model.second
		if !(:r in  keys(m))
			m[:r] = rterms(m[:vam], m[:data])
		end
		modtest.models[model.first] = m
	end
end

insert!(modtest, 
	:WInf => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	), 
	:W1 => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
		# 	"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
		# ] 
	),
	:WInf_bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)))
	),
	:W1_bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)))
	),
	:LLInf => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	), 
	:LL1 => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | LogLinear(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
		# 	"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
		# ] 
	),
	:LLInf_bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5)))
	),
	:LL1_bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA1(0.4) | LogLinear(0.001,2.5)))
	)

)

modtest.models[:WInf][:r]


function update!(modtest::ModelTest, key::Symbol)
	model = modtest.models[key]
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
	modtest.results[key]=result
end

function update!(modtest::ModelTest)
	for (key, _) in modtest.models
		update!(modtest, key)
	end
end

function test(modtest::ModelTest)
	for (key, result) in modtest.results
		@testset verbose = true "model $key" begin
			@testset "result $k" for k in [:lnL, :dlnL, :C, :dC]
				@test result[:resJL][k] ≈ result[:resR][k] atol=0.00000000000001
			end
		end
	end
end

function empty!(modtest::ModelTest; mode=:all)
	if mode ∈ [:all, :models]
		empty!(modtest.models)
	end
	if mode ∈ [:all, :results]	
		empty!(modtest.results)
	end
end


#empty!(modtest)
update!(modtest, :WInf)
update!(modtest)
test(modtest)

modtest.results[:WInf2][:resJL][:dC] .- modtest.results[:WInf2][:resR][:dC]
# string(modtest.models[:W1][:vam].formula)
# rterms(modtest.models[:W1][:vam], DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]))