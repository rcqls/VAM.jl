struct Sim
    model::Model
    stop_policy::StopPolicy
end

function init(sim::Sim, cache_size::Int)
    #// Almost everything in the 5 following lines are defined in model->init_computation_values() (but this last one initializes more than this 5 lines)
    sim.model.Vright=0
    sim.model.A=1
    sim.model.k=0
    for i in 0:(sim.model.nbPM + 1) 
        init(sim.model.models[i])
    end
    # size=cache_size_+1;cache_size=cache_size_;
    sim.model.idMod=0 #// Since no maintenance is possible!
    # //model->time=rep(0,size);
    # //model->type= rep(1,size);
    empty!(model.time)
    empty!(model.type)
    # (model->time).resize(size,0);
    # (model->type).resize(size,1);
end

function simulate(sim::Sim, nbsim::Int)::DataFrame
    init(sim, nbsim)

    first(sim.stop_policy)

    has_maintenance_policy(sim.model) || first(sim.model.maintenance_policy)

    while ok(sim, sim.stop_policy)
        #printf("k=%d\n",model->k);
        #To dynamically increase the size of simulation
        #resize();
        #printf("k2=%d\n",model->k);

        # //### modAV <- if(Type[k]<0) obj$vam.CM[[1]]$model else obj$vam.PM$models[[obj$data$Type[k]]]
        # //# Here, obj$model$k means k-1
        # //#print(c(obj$model$Vleft,obj$model$Vright))
        u = log(rand(1)[1])::Float64
        if sim.model.nb_paramsCov > 0
            u *= compute_covariates(sim) #;//set_current_system launched in R for simulation
        end
        timePM= 0.0
        #timeCM = sim.model->virtual_age_inverse(sim.model->family->inverse_cumulative_hazardRate(sim.model->family->cumulative_hazardRate(sim.model->virtual_age(sim.model->time[sim.model->k]))-u))
        timeCM = virtual_age_inverse(sim.model, inverse_cumulative_hazard_rate(sim.model.family, cumulative_hazard_rate(sim.model.family, virtual_age(sim.model,sim.model.time[sim.model.k]))-u))
       
        #//TODO: submodels
        idMod = 0
        # List timeAndTypePM;
        if has_maintenance_policy(sim.model)
            timeAndTypePM = update(sim.model.maintenance_policy, sim.model) # //# Peut-être ajout Vright comme argument de update
            # //timeAndTypePM = model->maintenance_policy->update(model->time[model->k]); //# Peut-être ajout Vright comme argument de update

            # //DEBUG[distrib type1]:
            # //NumericVector tmp0=timeAndTypePM["type"];
            # //int type0PM=tmp0[0];
            # //typeCptAV++;if(type0PM==1) type1CptAV++;

            # NumericVector tmp=timeAndTypePM["time"];
            # timePM=tmp[0];
            timePM = timeAndTypePM.time[1]
            #//DEBUG[distrib type1]: printf("sim: timePM(%d):%lf, timeCM=%lf\n",type0PM,timePM,timeCM);
            if timePM < timeCM && timePM < sim.model.time[model->k]
          		#print("Warning: PM ignored since next_time(=%lf)<current_time(=%lf) at rank %d.\n",timePM,model->time[model->k],model->k);
                print("warning")
            end
        end
        if !has_maintenance_policy(sim.model) || timeCM < timePM || timePM<sim.model.time[sim.model.k]
            sim.model.time[sim.model.k + 1]=timeCM
            sim.model.type[sim.model.k + 1]=-1
            idMod=0
        else
            sim.model.time[sim.model.k + 1]=timePM
            # NumericVector tmp2=timeAndTypePM["type"];
            # int typePM=tmp2[0];
            typePM = timeAndTypePM.type
            #//DEBUG[distrib type1]: typeCptAP++;if(typePM==1) type1CptAP++;printf("typePM=%d\n",typePM);
            sim.model.type[sim.model.k + 1]=typePM
            idMod=timeAndTypePM.type
        end
        #//printf("k=%d: cm=%lf,pm=%lf\n",model->k,timeCM,timePM);
        #//# used in the next update
        update_Vleft(sim.model, false,false)


        #//# update the next k, and save model in model too!
        update(sim.model.models[idMod],false,false)

    end
    #//DEBUG[distrib type1]: printf("cpt: %d/%d et %d/%d\n",type1CptAV,typeCptAV,type1CptAP,typeCptAP);

    return get_last_data(sim)
end

## Ok methods related to StopPolicy
ok(sim::Sim, sp::SizeGreaterThanStopPolicy)::Bool = sim.model.k < sp.size

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

ok(sim::Sim, sp::TimeGreaterThanStopPolicy)::Bool = sim.model.time[sim.model.k] < sp.time

function ok(sim::Sim, sp::AndStopPolicy)::Bool
    ans = false
    for policy in sp.policies 
        if ok(sim, policy)
            ans |= true; #every cond is tested because some init is done there!
        end
    end
    return ans
end

function ok(sim::Sim, sp::OrStopPolicy)::Bool
    ans = false
    for policy in sp.policies 
        if ok(sim, policy)
            ans &= false; #every cond is tested because some init is done there!
        end
    end
    return ans
end