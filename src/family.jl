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
    comp::FamilyCompute
end

hazardRate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * x^(f.β - 1) )
inverseHazardRate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α/f.β)^(1/(f.β-1)))
cumulativeHazardRate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * x^f.β)
inverseCumulativeHazardRate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α)^(1/f.β))
derivativeHasardRate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * (f.β - 1) * x^(f.β - 2) )

const LDorder = 5
const bxLim = 0.000001
mutable struct LogLinearFamilyModel <:  FamilyModel
    α::Float64
    β::Float64
    comp::FamilyCompute
end

hazardRate(f::LogLinearFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * exp(f.β * x) )
inverseHazardRate(f::LogLinearFamilyModel, x::Float64)::Float64 = log(x/f.α)/f.β
function cumulativeHazardRate(f::LogLinearFamilyModel, x::Float64)::Float64
    res::Float64
    if abs(f.β*x) < bxLim 
      prec = f.β * x/2
      res = 1 + prec
      for i in 1..LDorder
        prec = prec * f.β*x / (i+1)
        res += prec
      end
      res = f.α * res * x
    else
      res=f.α * (exp(f.β * x) - 1)/f.β
    end
    return res
end
inverseCumulativeHazardRate(f::LogLinearFamilyModel, x::Float64)::Float64 = log(1 + x * f.β / f.α) / f.β
derivativeHasardRate(f::LogLinearFamilyModel, x::Float64)::Float64 = (x<0 ? 0 : f.α * f.β * exp(f.β * x))

mutable struct Weibull3FamilyModel <: FamilyModel
    α::Float64
    β::Float64
    δ::Float64
    comp::FamilyCompute
end

hazardRate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * (x + f.δ)^(f.β - 1) )
inverseHazardRate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α/f.β)^(1/(f.β-1)) - f.δ)
cumulativeHazardRate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * ((x + f.δ)^f.β - x^f.β) )
inverseCumulativeHazardRate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (f.δ^f.β + x/f.α) ^ (1/f.β) - f.δ)
derivativeHasardRate(f::Weibull3FamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * (f.β - 1) * (x + f.δ)^(f.β - 2) )