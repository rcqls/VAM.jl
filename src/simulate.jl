struct Sim
    model::Model
    stop_policy::StopPolicy
end

function init(sim::Sim, cache_size::Int)
    #// Almost everything in the 5 following lines are defined in model->init_computation_values() (but this last one initializes more than this 5 lines)
    sim.model.Vright=0
    sim.model.A=1
    sim.model.k=0
    for i in 0..(model->nbPM + 1) 
        init(sim.model.models[i])
    end
    size=cache_size_+1;cache_size=cache_size_;
    model->idMod=0; #// Since no maintenance is possible!
    # //model->time=rep(0,size);
    # //model->type= rep(1,size);
    empty!(model->time)
    empty!(model->type)
    # (model->time).resize(size,0);
    # (model->type).resize(size,1);
end

function simulate(sim::Sim ,nbsim::Int)::DataFrame
    init(sim, nbsim)

    first(sim.stop_policy)

    has_maintenance_policy(sim.model) || first(sim.model.maintenance_policy)

    while(ok(sim.stop_policy))
        #printf("k=%d\n",model->k);
        #To dynamically increase the size of simulation
        #resize();
        #printf("k2=%d\n",model->k);

        # //### modAV <- if(Type[k]<0) obj$vam.CM[[1]]$model else obj$vam.PM$models[[obj$data$Type[k]]]
        # //# Here, obj$model$k means k-1
        # //#print(c(obj$model$Vleft,obj$model$Vright))
        u=log(rand(1)[1])::Float64
        if(sim.model.nb_paramsCov > 0) u *= compute_covariates(sim) #;//set_current_system launched in R for simulation
        timePM= 0.0
        #timeCM = sim.model->virtual_age_inverse(sim.model->family->inverse_cumulative_hazardRate(sim.model->family->cumulative_hazardRate(sim.model->virtual_age(sim.model->time[sim.model->k]))-u))
        timeCM = virtual_age_inverse(sim.model, inverse_cumulative_hazardRate(sim.model->family, cumulative_hazardRate(sim.model->family, virtual_age(sim.model,sim.model->time[sim.model->k]))-u))
       
        #//TODO: submodels
        idMod = 0
        # List timeAndTypePM;
        if(has_maintenance_policy(sim.model))
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
            if(timePM < timeCM && timePM < sim.model.time[model->k])
          		#print("Warning: PM ignored since next_time(=%lf)<current_time(=%lf) at rank %d.\n",timePM,model->time[model->k],model->k);
                print("warning")
            end
        end
        if(!has_maintenance_policy(sim.model) || timeCM < timePM || timePM<sim.model.time[sim.model.k])
            sim.model.time[sim.model.k + 1]=timeCM
            sim.model.type[sim.model.k + 1]=-1
            idMod=0
        else
            sim.model.time[sim.model.k + 1]=timePM
            # NumericVector tmp2=timeAndTypePM["type"];
            # int typePM=tmp2[0];
            typePM = timeAndTypePM.type[1]
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

# function add_stop_policy(List policy) {
#     stop_policy=newStopPolicy(this,policy);
# }

# void SimVam::init(int cache_size_) {
#     int i;
#     // Almost everything in the 5 following lines are defined in model->init_computation_values() (but this last one initializes more than this 5 lines)
#     model->Vright=0;
#     model->A=1;
#     model->k=0;
#     for(i=0;i<model->nbPM + 1;i++) model->models->at(i)->init();

#     size=cache_size_+1;cache_size=cache_size_;
#     model->idMod=0; // Since no maintenance is possible!
#     //model->time=rep(0,size);
#     //model->type= rep(1,size);
#     (model->time).clear();
#     (model->type).clear();
#     (model->time).resize(size,0);
#     (model->type).resize(size,1);
# }
# public:

#     SimVam(List model_) {
#         model=new VamModel(model_);
#     };

#     ~SimVam() {
#         delete model;
#     };

#     //This is used inside simulate to fetch last generated data
#     DataFrame get_last_data();

#     //The 2 following functions are now similarly used in plot.R than mle.vam and model.vam
#     //This is just a delagation to model 
#     void set_data(List data_) {
#         model->set_data(data_);
#     }

#     DataFrame get_selected_data(int i) {
#         return model->get_selected_data(i);
#     }

#     DataFrame simulate(int nbsim);

#     VamModel* get_model() {
#         return model;
#     }

#     NumericVector get_params() {
#         return model->get_params();
#     }

#     void set_params(NumericVector pars) {
#         model->set_params(pars);
#     }

#     //delegate from model cache!
#     List get_virtual_age_infos(double by,double from, double to) {
#         return model->get_virtual_age_infos(by,from,to);
#     }

#     void add_stop_policy(List policy);

#     int cache_size,size;

#     //Covariates related

#     void select_current_system(int i) {
#         model->select_current_system(i,false);
#     }

# 	double compute_covariates() {
#         return exp(-model->compute_covariates());
#     }; 

#     void set_covariates(List covariates_) {
#         model->set_covariates(covariates_);
#     }

# private:
#     VamModel* model;

#     StopPolicy* stop_policy;

#     void init(int cache_size_);

#     void resize();

# };

#endif //RCPP_SIM_VAM_H
