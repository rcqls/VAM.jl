function simulate(model::Model, stop::Union{Nothing, Real, Vector{Any}}; system::Int=1)::DataFrame
    sim = simulator(model, stop)
    return simulate(sim, system = system)
end

mutable struct Simulator
    model::Model
    stop_policy::Union{Nothing,Expr}
end

function simulator(model::Model, stop::Union{Nothing, Real, Vector{Any}})::Simulator
    sim = Simulator(model, nothing)
    add_stop_policy!(sim, stop)
    init!(sim.model)
    return sim
end

sim(model::Model, stop::Union{Nothing, Real, Vector{Any}})::Simulator = simulator(model, stop)

function init!(sim::Simulator)
    #// Almost everything in the 5 following lines are defined in model->init_computation_values() (but this last one initializes more than this 5 lines)
    sim.model.Vright=0
    sim.model.A=1
    sim.model.k=1
    for mm in sim.model.models
        init!(mm)
    end
    # size=cache_size_+1;cache_size=cache_size_;
    sim.model.id_mod=0 #// Since no maintenance is possible!
    sim.model.time = [0.0]
    sim.model.type = [-1]
end

function simulate(sim::Simulator, stop::Union{Nothing, Real, Vector{Any}}; system::Int=1)::DataFrame
    add_stop_policy!(sim, stop)
    if has_maintenance_policy(sim.model)
        first(sim.model.maintenance_policy)
    end
    data = DataFrame()
    for syst in 1:system
        init!(sim)
        run = true
        while run
            u = log(rand(1)[1])::Float64
            if sim.model.nb_params_cov > 0
            #   u *= compute_covariates(sim) #;//set_current_system launched in R for simulation
            end
            timeCM = virtual_age_inverse(sim.model, inverse_cumulative_hazard_rate(sim.model.family, cumulative_hazard_rate(sim.model.family, virtual_age(sim.model,sim.model.time[sim.model.k]))-u))
            ##println("timeCM=$timeCM")
            #   TODO: submodels
            id_mod = 0
        #     # List timeAndTypePM;
            if has_maintenance_policy(sim.model)
                timePM, typePM = update(sim.model.maintenance_policy, sim.model) # //# Peut-être ajout Vright comme argument de update
                if timePM < timeCM && timePM < sim.model.time[sim.model.k]
                    #print("Warning: PM ignored since next_time(=%lf)<current_time(=%lf) at rank %d.\n",timePM,model->time[model->k],model->k);
                    print("warning")
                end
            end
            if !has_maintenance_policy(sim.model) || timeCM < timePM || timePM < sim.model.time[sim.model.k]
                push!(sim.model.time,timeCM)
                push!(sim.model.type, -1)
                id_mod=0
            else
                push!(sim.model.time, timePM)
                #//DEBUG[distrib type1]: typeCptAP++;if(typePM==1) type1CptAP++;printf("typePM=%d\n",typePM);
                push!(sim.model.type, typePM)
                id_mod=typePM
            end
        #     #//printf("k=%d: cm=%lf,pm=%lf\n",model->k,timeCM,timePM);
        #     #//# used in the next update
            update_Vleft!(sim.model) #, false,false)


        #     #//# update the next k, and save model in model too!
            update_maintenance!(sim.model, id_mod) #false,false)
            run = ok(sim)
            ## TODO work on stop later
        end
        data = vcat(data,DataFrame(system=syst, time=sim.model.time, type=sim.model.type))
    end
    if system == 1
        data = data[:,[:time, :type]]
    end
    df = data[2:size(data,1),:]
    if system > 1
        rename!(df, vcat(["System"], sim.model.varnames) )
    else
        rename!(df, sim.model.varnames)
    end
    df
end

simulate(sim::Simulator; system::Int=1) = simulate(sim, nothing, system=system)

function ok(sim::Simulator)::Bool
    s = length(sim.model.time) - 1 # 1st is 0 time to be removed when returned
    t = sim.model.time[sim.model.k]
    eval(:(s=$s))
    eval(:(t=$t))
    eval(sim.stop_policy)
end

function add_stop_policy!(sim::Simulator, stop::Union{Nothing, Real,Vector{Any}})
    if isa(stop,Int)
        sim.stop_policy = Expr(:call, :<=, :s,  stop)
    elseif isnothing(stop)
        if isnothing(sim.stop_policy)
            sim.stop_policy = Expr(:call, :<=, :s,  100)
        end
    else  
        sim.stop_policy =  formula_translate(Expr(:call, stop...))
    end
end