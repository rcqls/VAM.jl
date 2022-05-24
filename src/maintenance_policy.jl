abstract type AbstractMaintenancePolicy end
abstract type MaintenancePolicyWithExternalModel <: AbstractMaintenancePolicy end

struct PeriodicMaintenancePolicy <: AbstractMaintenancePolicy
    from::Float64
    by::Float64
    prob::Vector{Float64}
end

PeriodicMaintenancePolicy(by::Real, prob::Vector{Float64}) = PeriodicMaintenancePolicy(0.0, by, prob)
PeriodicMaintenancePolicy(from::Real, by::Real) = PeriodicMaintenancePolicy(from, by, [1])
PeriodicMaintenancePolicy(by::Real) = PeriodicMaintenancePolicy(0, by)

type_size(mp::PeriodicMaintenancePolicy)::Int = length(mp.prob)

function first(mp::PeriodicMaintenancePolicy); end

function update(mp::PeriodicMaintenancePolicy, model::AbstractModel)::NamedTuple{(:time, :type), Tuple{Float64, Int64}}
    current=model.time[model.k]
	time = mp.from + (floor((current - mp.from)/mp.by) + 1) * mp.by
 
	r=rand(1)[1]
    t = 0
	for p in prob[1:end-1]
        if r < p
            break
        else 
            t += 1
            r -= p
        end
    end
	# //printf("from=%d,t=%d, %lf\n",get_from_type(),t,prob[0]);
	return (time=time, type=t)
end


struct AtTimesMaintenancePolicy <:  AbstractMaintenancePolicy
    times::Vector{Float64}
    i::Int
    k::Int
    cycle::Bool
    differentTypeIfCM::Bool
end

type_size(mp::AtTimesMaintenancePolicy)::Int = mp.differentTypeIfCM ? 2 : 1

struct AtIntensityMaintenancePolicy <: MaintenancePolicyWithExternalModel
    level::Float64
    #external_model::AbstractModel
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
    from_type::Vector{Int}
end

function MaintenancePolicyList(policies::Vector{AbstractMaintenancePolicy})
    from_type = [0]
    from = 0
    for policy in policies[1:end - 1]
        from += type_size(policy)
        push!(from_type, from)
    end
    MaintenancePolicyList(policies, from_type)
end

function type_size(mp::MaintenancePolicyList)::Int
    s = 0
    for policy in mp.policies
        s += type_size(policy)
    end
    return s;
end