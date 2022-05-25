abstract type StopPolicy end

init!(sp::StopPolicy) = nothing

struct SizeGreaterThanStopPolicy <: StopPolicy
    size::Int
end
# //M is for CM or PM, type=-1 or 1,2,...

struct SizeOfTypeGreaterThanStopPolicy <:StopPolicy
    size::Int
    type::Int
    count::Int
end

### TODO LATER
# struct TimeGreaterThanCensorshipStopPolicy <: StopPolicy
#     to_init::Bool
#     time::Float64
#     #expr::Language
#     #env::Environment
# end

# function ok(sim::Sim, sp::TimeGreaterThanCensorshipStopPolicy)::Bool
#     ok = sim.model.time[sim.model.k] < sp.time
#     if !ok
#         #//update result of sim!
#         sim.model.time[sim.model.k]=sp.time
#         si.model.type[sim.model.k]=0
#     end
#     return ok
# end

# function TimeGreaterThanCensorshipStopPolicy::first() {
#     //printf("time=%lf\n",time);
#      if(to_init) {
#        time=as<double>(Rf_eval(expr,env));
#      };
# }

struct TimeGreaterThanStopPolicy <: StopPolicy
    time::Float64
end
struct AndStopPolicy <: StopPolicy
    policies::Vector{StopPolicy}
end

struct OrStopPolicy <: StopPolicy
    policies::Vector{StopPolicy}
end
 