abstract type FamilyModel end

## TODO: FamilyCompute NOT USED from now, used in the non yet implemented  
# hazardRate_param_derivative 
# cumulative_hazardRate_param_derivative 
# hazardRate_derivative_param_derivative 
# hazardRate_2derivative 
# hazardRate_param_2derivative 
# cumulative_hazardRate_param_2derivative

# mutable struct FamilyCompute
#     dHR::Vector{Float64}
#     dHL::Vector{Float64}
#     dhR::Vector{Float64}
#     dhL::Vector{Float64}
#     d2HR::Matrix{Float64}
#     d2HL::Matrix{Float64}
#     d2h::Matrix{Float64}
#     dhd::Matrix{Float64}
#     nb_params::Int
# end
mutable struct WeibullFamilyModel <: FamilyModel
    α::Float64
    β::Float64
    covariates::Union{Nothing, Vector{Any}}
    # comp::FamilyCompute
end
WeibullFamilyModel(α::Float64, β::Float64) = WeibullFamilyModel(α, β, nothing)
params(fm::WeibullFamilyModel)::Vector{Float64} = [fm.α,fm.β]
params!(fm::WeibullFamilyModel, p::Vector{Float64}) = begin; fm.α,fm.β = p; nothing; end
nb_params(fm::WeibullFamilyModel) = 2

hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * x^(f.β - 1) )
inverse_hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α/f.β)^(1/(f.β-1)))
cumulative_hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * x^f.β)
inverse_cumulative_hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α)^(1/f.β))
derivative_hasard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * (f.β - 1) * x^(f.β - 2) )

const LDorder = 5
const bxLim = 0.000001
mutable struct LogLinearFamilyModel <:  FamilyModel
    α::Float64
    β::Float64
    covariates::Union{Nothing, Vector{Any}}
    # comp::FamilyCompute
end
params(fm::LogLinearFamilyModel)::Vector{Float64} = [fm.α,fm.β]
params!(fm::LogLinearFamilyModel, p::Vector{Float64}) = begin;fm.α,fm.β = p; nothing; end
nb_params(fm::LogLinearFamilyModel) = 2

hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * exp(f.β * x) )
inverse_hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64 = log(x/f.α)/f.β
function cumulative_hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64
    res::Float64
    if abs(f.β*x) < bxLim 
      prec = f.β * x/2
      res = 1 + prec
      for i in 1:LDorder
        prec = prec * f.β*x / (i+1)
        res += prec
      end
      res = f.α * res * x
    else
      res=f.α * (exp(f.β * x) - 1)/f.β
    end
    return res
end
inverse_cumulative_hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64 = log(1 + x * f.β / f.α) / f.β
derivative_hasard_rate(f::LogLinearFamilyModel, x::Float64)::Float64 = (x<0 ? 0 : f.α * f.β * exp(f.β * x))

mutable struct Weibull3FamilyModel <: FamilyModel
    α::Float64
    β::Float64
    δ::Float64
    covariates::Union{Nothing, Vector{Any}}
    # comp::FamilyCompute
end
WeibullFamilyModel(α::Float64, β::Float64, δ::Float64) = Weibull3FamilyModel(α, β, δ, nothing)
params(fm::Weibull3FamilyModel)::Vector{Float64} = [fm.α,fm.β,fm.δ]
params!(fm::Weibull3FamilyModel, p::Vector{Float64}) = begin; fm.α,fm.β,fm.δ = p; nothing; end
nb_params(fm::Weibull3FamilyModel) = 3

hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * (x + f.δ)^(f.β - 1) )
inverse_hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α/f.β)^(1/(f.β-1)) - f.δ)
cumulative_hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * ((x + f.δ)^f.β - x^f.β) )
inverse_cumulative_hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (f.δ^f.β + x/f.α) ^ (1/f.β) - f.δ)
derivative_hasard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * (f.β - 1) * (x + f.δ)^(f.β - 2) )