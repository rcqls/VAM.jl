abstract type AbstractMaintenancePolicy end
abstract type MaintenancePolicyWithExternalModel <: AbstractMaintenancePolicy end

struct PeriodicMaintenancePolicy <: AbstractMaintenancePolicy
    from::Vector{Float64}
    by::Vector{Float64}
    prob::Vector{Float64}
end
    
type_size(mp::PeriodicMaintenancePolicy)::Int = length(mp.prob)
function first(mp::PeriodicMaintenancePolicy); end

struct AtTimesMaintenancePolicy <:  AbstractMaintenancePolicy
    times::Vector{Float64}
    i::Int
    k::Int
    cycle::Bool
    differentTypeIfCM::Bool
end

type_size(mp::AtTimesMaintenancePolicy)::Int = mp.differentTypeIfCM ? 2 : 1

struct AtIntensityMaintenancePolicy <: MaintenancePolicyWithExternalModel
    level::Vector{Float64}
    external_model::AbstractModel
end

type_size(mp::AtIntensityMaintenancePolicy)::Int = 1


struct AtVirtualAgeMaintenancePolicy <:  MaintenancePolicyWithExternalModel
    level::Vector{Float64}
    external_model::AbstractModel
end

type_size(mp::AtVirtualAgeMaintenancePolicy)::Int = 1


struct AtFailureProbabilityMaintenancePolicy <: AbstractMaintenancePolicy
    level::Vector{Float64}
    external_model::AbstractModel
end

type_size(mp::AtFailureProbabilityMaintenancePolicy)::Int = 1
 

struct MaintenancePolicyList <: AbstractMaintenancePolicy
    policies::Vector{AbstractMaintenancePolicy}
end

function type_size(mp::MaintenancePolicyList)::Int
    s = 0
    for policy in mp.policies
        s += type_size(policy)
    end
    return s;
end