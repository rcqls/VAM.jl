mutable struct Sim
    model::Model
    stop_policy::Union{Nothing,Expr}
end

function init!(sim::Sim)
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

function simulate(sim::Sim, stop::Union{Nothing, Real, Vector{Any}}; system::Int=1)::DataFrame
    add_stop_policy!(sim, stop)
    has_maintenance_policy(sim.model) || first(sim.model.maintenance_policy)
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
            update!(sim.model.models[id_mod + 1], sim.model) #false,false)
            save_id_mod(sim.model, id_mod)
            run = ok(sim)
            ## TODO work on stop later
        end
        data = vcat(data,DataFrame(system=syst, time=sim.model.time, type=sim.model.type))
    end
    if system ==1
        data = data[:,[:time, :type]]
    end
    data
end

simulate(sim::Sim; system::Int=1) = simulate(sim, nothing, system=system)

function ok(sim::Sim)::Bool
    s = length(sim.model.time)
    t = sim.model.time[sim.model.k]
    eval(:(__size__=$s;__time__=$t))
    eval(sim.stop_policy)
end

function add_stop_policy!(sim::Sim, stop::Union{Nothing, Real,Vector{Any}})
    if isa(stop,Int)
        sim.stop_policy = Expr(:call, :<=, :__size__,  stop)
    elseif isnothing(stop)
        if isnothing(sim.stop_policy)
            sim.stop_policy = Expr(:call, :<=, :__size__,  100)
        end
    else  
        sim.stop_policy = formula_translate(Expr(:call, stop...))
    end
end