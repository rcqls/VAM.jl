abstract type StopPolicy end

function first(sp::StopPolicy); end

struct SizeGreaterThanStopPolicy <: StopPolicy
    size::Int
end

ok(sim::Sim, sp::SizeGreaterThanStopPolicy)::Bool = sim.model.k < sp.size


# //M is for CM or PM, type=-1 or 1,2,...

struct SizeOfTypeGreaterThanStopPolicy <:StopPolicy
    size::Int
    type::Int
    count::Int
end

function ok(sim::Sim, sp::SizeOfTypeGreaterThanStopPolicy)::Bool
    # //incr counter
    # //printf("k=%d,t=%lf,ty=%d ->",mod->k,mod->time[mod->k],mod->type[mod->k]);
    if sim.model.k>0 && sim.model.type[sim.model.k] == sp.type
        sp.count += 1
        #//printf("type=%d,count=%d",type,count);
    end
    #//printf("\n");
    return sp.count < sp.size
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

ok(sim::Sim, sp::TimeGreaterThanStopPolicy)::Bool = sim.model.time[sim.model.k] < sp.time

struct AndStopPolicy <: StopPolicy
    policies::Vector{StopPolicy}
end

struct OrStopPolicy <: StopPolicy
    policies::Vector{StopPolicy}
end

function ok(sim::Sim, sp::AndStopPolicy)::Bool
    ans = false
    for policy in sp.policies 
        if ok(sim, policy)
            ans |= true; #every cond is tested because some init is done there!
        end
    end
    return ans
end

function ok(sim::Sim, sp::ORStopPolicy)::Bool
    ans = false
    for policy in sp.policies 
        if ok(sim, policy)
            ans &= false; #every cond is tested because some init is done there!
        end
    end
    return ans
end


# bool AndStopPolicy::ok() {
#     bool ans=false;
#     for(
#         std::vector<StopPolicy*>::iterator it=policies.begin();
#         it != policies.end();
#         ++it
#     ) {
#         if((*it)->ok()) ans |= true; //every cond is tested because some init is done there!
#     }
#     return ans;
# }

# bool OrStopPolicy::ok() {
#     bool ans=true;
#     for(
#         std::vector<StopPolicy*>::iterator it=policies.begin();
#         it != policies.end();
#         ++it
#     ) {
#         if(!(*it)->ok()) ans &= false; //every cond is tested because some init is done there!
#     }
#     return ans;
# }


# // Exactly the same first method for the following classes
# void AndStopPolicy::first() {
#     for(
#         std::vector<StopPolicy*>::iterator it=policies.begin();
#         it != policies.end();
#         ++it
#     ) {
#         (*it)->first();
#     }
# }

# void OrStopPolicy::first() {
#     for(
#         std::vector<StopPolicy*>::iterator it=policies.begin();
#         it != policies.end();
#         ++it
#     ) {
#         (*it)->first();
#     }
# }