mutable struct MLE
    model::Model
    left_censors::Vector{Int} #CAREFUL: this is a vector of indices!
    left_censor::Int #left_censor for current system

    comp::Compute
    MLE() = new()
end
function MLE(model::Model, data::DataFrame)::MLE
    mle = MLE()
    mle.model = model
    init!(mle.model)
    mle.comp = Compute(mle.model)
    data!(mle.model, data)
    left_censors!(mle, Int[])
    return mle
end

mle(model::Model, data::DataFrame)::MLE = MLE(model, data)

## TODO: deal with left_censors
function left_censors!(m::MLE, left_censors::Vector{Int})
      m.left_censors = left_censors
end

function select_left_censor(mle::MLE, i::Int)
    if length(mle.left_censors) > 0 
        mle.left_censor=left_censors[i]
    end
end

function contrast(mle::MLE, param::Vector{Float64}; alpha_fixed::Bool=true)::Float64
    res = 0
    alpha = param[1] #;//save current value of alpha

    param[1] = 1 #//Rmk: alpha replaces param[0] => a bit weird!

    init!(mle.comp)

    params!(mle.model,param);
    # //printf("System %d\n",1);
    select_data(mle.model, 1)
    # if(model.nb_params_cov > 0) model.select_current_system(0,true);
    select_left_censor(mle, 1)
    contrast_current(mle)
    # //only if multi-system
    if mle.model.nb_system > 1
        for i in 2:mle.model.nb_system 
            #//printf("System %d\n",i+1);
            select_data(mle.model, i)
            if model.nb_params_cov > 0 
                select_current_system(mle.model, i, true)
            end
            select_left_censor(mle, i)
            contrast_for_current_system()
        end
    end

    # //DEBUG: printf("alpha=%lf,S0=%lf,S1=%lf,S2=%lf,S3=%lf,S4=%lf\n",alpha,S0,S1,S2,S3,S4);
    # //printf("params=(%lf,%lf)\n",model.params_cov[0],model.params_cov[1]);
    # // log-likelihood (with constant +S0*(log(S0)-1))
    if !alpha_fixed
        param[1] = mle.comp.S0 / mle.comp.S1
        res = -log(mle.comp.S1) * mle.comp.S0 + mle.comp.S2 +  mle.comp.S0 * (log( mle.comp.S0) - 1) +  mle.comp.S3;
        params!(mle.model, param) #//also memorize the current value for alpha which is not 1 in fact
    else
        param[1] = alpha
        res = log(alpha) * mle.comp.S0 + mle.comp.S2 - alpha * mle.comp.S1 + mle.comp.S3
        params!(mle.model, param) #//also memorize the current value for alpha which is not 1 in fact
    end
    if mle.model.nb_params_cov > 0 
        res += mle.comp.S4
    end

    param[1] = alpha #//LD:changed for bayesian
    return res
end

function contrast_current(mle::MLE)
    init_mle(mle)
    n = length(mle.model.time)
    while mle.model.k < n
        # //printf("  Time=%f, Type=%d\n",model.time[model.k+1],model.type[model.k+1]);
        contrast_update_current(mle)
        #// previous model for the next step
        type = mle.model.type[mle.model.k + 1]
        if type < 0 
            type = 0
        end
        #//model.indMode = (type < 0 ? 0 : type);
        update!(mle.model.models[1 + type], mle.model)
    end
    contrast_update_S(mle)
end

function contrast_update_current(mle::MLE; deriv::Bool=false)
    update_Vleft!(mle.model,with_gradient=deriv,with_hessian=deriv)
    mle.model.hVleft = hazard_rate(mle.model.family, mle.model.Vleft)
    mle.model.indType = (mle.model.type[mle.model.k + 1] < 0 ? 1.0 : 0.0)
    if mle.model.k >= mle.left_censor 
        mle.model.comp.S1 += cumulative_hazard_rate(mle.model.family, mle.model.Vleft) - cumulative_hazard_rate(mle.model.family, mle.model.Vright)
    end
    mle.model.comp.S2 += log(mle.model.hVleft) * mle.model.indType
    mle.model.comp.S3 += log(mle.model.A) * mle.model.indType
end

function contrast_update_S(mle::MLE)
    #//model updated for current system: S1,S2,S0,dS1,dS2
    tmp = mle.model.comp.S1
    mle.comp.S2 += mle.model.comp.S2
    mle.comp.S0 += mle.model.comp.S0
    mle.comp.S3 += mle.model.comp.S3
    if mle.model.nb_params_cov > 0
        compute_covariates(mle.model) #//initialize model.sum_cov
        tmp *= exp(mle.model.sum_cov)
        mle.comp.S4 += mle.model.comp.S0 * mle.model.sum_cov
        #//printf("(S0=%lf) * (sum_cov=%lf) = (S4 =%lf)\n",model.S0, model.sum_cov,model.S0 * model.sum_cov);
    end
    mle.comp.S1 += tmp
    #//printf("Conclusion : S1=%f, S2=%f, S0=%f, S4=%f\n",model.S1,model.S2,model.S0,model.S4);
end



#     void gradient_for_current_system() {
#         int i,ii;
#     	init_mle_vam_for_current_system(true,false);
#     	int n=(model.time).size() - 1;
#     	while(model.k < n) {
#     		gradient_update_for_current_system();
#     		int type=model.type[model.k + 1 ];
# 			if(type < 0) type=0;
# 			//model.indMode = (type < 0 ? 0 : type);
# 			model.models.at(type).update(true,false);
#     	}
#         contrast_S_update();
#         //precomputation of covariate term to multiply (in fact just exp)
#         for(i=0;i<model.nb_paramsFamily-1;i++) {
#             gradient_dS_family_update(i);
#         }
#         for(ii=0;ii<model.nb_params_maintenance;ii++,i++) {
#             gradient_dS_maintenance_update(i,ii);
#         }
#         for(ii=0;ii<model.nb_params_cov;ii++,i++) {
#             gradient_dS_covariate_update(i,ii);
#         }
#     }

#     NumericVector gradient(NumericVector param, bool alpha_fixed=false) {
#         NumericVector res(model.nb_paramsFamily+model.nb_params_maintenance+model.nb_params_cov);
#         double alpha=param[0];//save current value of alpha
#         int i,ii;

#         param[0]=1;
#         init_mle_vam(true,false);
#         model.set_params(param);
#         model.select_data(0);
#         if(model.nb_params_cov > 0) model.select_current_system(0,true);
#         select_left_censor(0);
#         gradient_for_current_system();
#         //only if multi-system
#         for(i=1;i<model.nb_system;i++) {
#             model.select_data(i);
#             if(model.nb_params_cov > 0) model.select_current_system(i,true);
#             select_left_censor(i);
#             gradient_for_current_system();
#         }
#         //compute gradient
#         param[0] = (alpha_fixed ? alpha : S0 / S1);

#         model.set_params(param);//also memorize the current value for alpha which is not 1 in fact

#         res[0] = (alpha_fixed ? S0/alpha-S1 : 0);
#         //
#         for(i=0;i<model.nb_paramsFamily-1;i++) {
#             res[i+1] = -dS1[i]*param[0] + dS2[i];
#         }
#         for(ii=0;ii<model.nb_params_maintenance;ii++,i++) {
#             res[i+1] = -dS1[i]*param[0] + dS2[i]+dS3[ii];
#         }
#         for(ii=0;ii<model.nb_params_cov;ii++,i++) {
#             res[i+1] = -dS1[i]*param[0] + dS4[ii] ;
#         }

#         param[0]=alpha;//LD:changed for bayesian
#         return res;
#     }

#     void hessian_for_current_system() {
#         int j,i,ii,k,kk;
#         init_mle_vam_for_current_system(true,true);
#         int n=(model.time).size() - 1;
#         while(model.k < n) {
#             hessian_update_for_current_system();
#             int type=model.type[model.k + 1 ];
#             if(type < 0) type=0;
#             //model.indMode = (type < 0 ? 0 : type);
#             model.models.at(type).update(true,true);
#         }
#         contrast_S_update();
#         for(i=0;i<model.nb_paramsFamily-1;i++) {
#             gradient_dS_family_update(i);
#             for(j=0;j<=i;j++) {
#                 //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                 k=i*(i+1)/2+j;
#                 d2S1[k] += model.d2S1[k] * (model.nb_params_cov > 0 ? exp(model.sum_cov) : 1.0);
#                 d2S2[k] += model.d2S2[k];
#             }
#         }
#         for(ii=0;ii<model.nb_params_maintenance;ii++,i++) {
#             gradient_dS_maintenance_update(i,ii);
#             for(j=0;j<=ii;j++) {
#                  //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                 k=i*(i+1)/2+j;
#                 d2S1[k] += model.d2S1[k] * (model.nb_params_cov > 0 ? exp(model.sum_cov) : 1.0);
#                 d2S2[k] += model.d2S2[k];
#                 //ii and j(<=ii) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                 kk=ii*(ii+1)/2+j;
#                 d2S3[kk] += model.d2S3[kk];
#             }
#             for(j=ii+1;j<=i;j++) {
#                 //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                 k=i*(i+1)/2+j;
#                 d2S1[k] += model.d2S1[k] * (model.nb_params_cov > 0 ? exp(model.sum_cov) : 1.0);
#                 d2S2[k] += model.d2S2[k];
#             }
#         }
#         for(ii=0;ii<model.nb_params_cov;ii++,i++) {
#             gradient_dS_covariate_update(i,ii);
#             for(j=0;j<model.nb_paramsFamily-1 + model.nb_params_maintenance;j++) {
#               k=i*(i+1)/2+j;
#               d2S1[k] += model.get_covariate(ii)* exp(model.sum_cov)*model.dS1[j];
#             }
#             for(j=model.nb_paramsFamily-1 + model.nb_params_maintenance;j<=i;j++){
#               k=i*(i+1)/2+j;
#               d2S1[k] += model.get_covariate(ii) * model.get_covariate(j - model.nb_paramsFamily+1 - model.nb_params_maintenance) * exp(model.sum_cov)*model.S1;
#             }
#         }
#     }

#     NumericMatrix hessian(NumericVector param, bool alpha_fixed=false) {
#         int j;
#         NumericMatrix res(model.nb_paramsFamily+model.nb_params_maintenance+model.nb_params_cov,model.nb_paramsFamily+model.nb_params_maintenance+model.nb_params_cov);
#         double alpha=param[0];//save current value of alpha

#         param[0]=1;
#         init_mle_vam(true,true);
#         model.set_params(param);
#         model.select_data(0);
#         if(model.nb_params_cov > 0) model.select_current_system(0,true);
#         select_left_censor(0);
#         hessian_for_current_system();

#         //only if multi-system
#         for(int i=1;i<model.nb_system;i++) {
#             model.select_data(i);
#             if(model.nb_params_cov > 0) model.select_current_system(i,true);
#             select_left_censor(i);
#             hessian_for_current_system();
#         }

#         //compute hessian
#         if(!alpha_fixed) {
#             res(0,0) = 0;
#             for(int i=0;i<model.nb_paramsFamily-1;i++) {
#                 res(0,i+1) = 0;
#                 res(i+1,0) = 0;
#                 res(i+1,i+1) = pow(dS1[i],2)/pow(S1,2) * S0-d2S1[i*(i+1)/2+i]/S1 * S0 + d2S2[i*(i+1)/2+i];
#                 for(j=0;j<i;j++) {
#                     //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                     res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S0-d2S1[i*(i+1)/2+j]/S1 * S0 + d2S2[i*(i+1)/2+j];
#                     res(j+1,i+1) = res(i+1,j+1);
#                 }
#             }
#             for(int i=(model.nb_paramsFamily-1);i<(model.nb_params_maintenance+model.nb_paramsFamily-1);i++) {
#                 res(0,i+1) = 0;
#                 res(i+1,0) = 0;
#                 res(i+1,i+1) = pow(dS1[i],2)/pow(S1,2) * S0-d2S1[i*(i+1)/2+i]/S1 * S0 + d2S2[i*(i+1)/2+i] + d2S3[(i-(model.nb_paramsFamily-1))*(i-(model.nb_paramsFamily-1)+1)/2+i-(model.nb_paramsFamily-1)];
#                 for(j=0;j<(model.nb_paramsFamily-1);j++) {
#                     //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                     res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S0-d2S1[i*(i+1)/2+j]/S1 * S0 + d2S2[i*(i+1)/2+j];
#                     res(j+1,i+1) = res(i+1,j+1);
#                 }
#                 for(j=(model.nb_paramsFamily-1);j<i;j++) {
#                     //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                     res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S0-d2S1[i*(i+1)/2+j]/S1 * S0 + d2S2[i*(i+1)/2+j]+ d2S3[(i-(model.nb_paramsFamily-1))*(i-(model.nb_paramsFamily-1)+1)/2+j-(model.nb_paramsFamily-1)];
#                     res(j+1,i+1) = res(i+1,j+1);
#                 }
#             }
#             for(int i=(model.nb_params_maintenance+model.nb_paramsFamily-1);i<(model.nb_params_maintenance+model.nb_paramsFamily+model.nb_params_cov-1);i++) {
#                res(0,i+1) = 0;
#                res(i+1,0) = 0;
#                res(i+1,i+1) = pow(dS1[i],2)/pow(S1,2) * S0-d2S1[i*(i+1)/2+i]/S1 * S0;
#                for(j=0;j<i;j++) {
#                  res(i+1,j+1) = dS1[i]*dS1[j]/pow(S1,2) * S0-d2S1[i*(i+1)/2+j]/S1 * S0;
#                  res(j+1,i+1) = res(i+1,j+1);               }
#             }
#             param[0]=S0/S1;
#             model.set_params(param);//also memorize the current value for alpha which is not 1 in fact
#         } else {

#             res(0,0) = -S0/pow(alpha,2);
#             for(int i=0;i<(model.nb_paramsFamily-1);i++) {
#                 res(0,i+1) = -dS1[i];
#                 res(i+1,0) = -dS1[i];
#                 res(i+1,i+1) = d2S2[i*(i+1)/2+i]-alpha*d2S1[i*(i+1)/2+i];
#                 for(j=0;j<i;j++) {
#                     //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                     res(i+1,j+1) = d2S2[i*(i+1)/2+j]-alpha*d2S1[i*(i+1)/2+j];
#                     res(j+1,i+1) = res(i+1,j+1);
#                 }
#             }
#             for(int i=(model.nb_paramsFamily-1);i<(model.nb_params_maintenance+model.nb_paramsFamily-1);i++) {
#                 res(0,i+1) = -dS1[i];
#                 res(i+1,0) = -dS1[i];
#                 res(i+1,i+1) = d2S2[i*(i+1)/2+i]-alpha*d2S1[i*(i+1)/2+i]+d2S3[(i-(model.nb_paramsFamily-1))*(i- (model.nb_paramsFamily-1)+1)/2+i-(model.nb_paramsFamily-1)];
#                 for(j=0;j<(model.nb_paramsFamily-1);j++) {
#                     //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                     res(i+1,j+1) = d2S2[i*(i+1)/2+j]-alpha*d2S1[i*(i+1)/2+j];
#                     res(j+1,i+1) = res(i+1,j+1);
#                 }
#                 for(j=(model.nb_paramsFamily-1);j<i;j++) {
#                     //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                     res(i+1,j+1) = d2S2[i*(i+1)/2+j]-alpha*d2S1[i*(i+1)/2+j]+d2S3[(i-(model.nb_paramsFamily-1))*(i- (model.nb_paramsFamily-1)+1)/2+j-(model.nb_paramsFamily-1)];
#                     res(j+1,i+1) = res(i+1,j+1);
#                 }
#             }
#             for(int i=(model.nb_params_maintenance+model.nb_paramsFamily-1);i<(model.nb_params_maintenance+model.nb_paramsFamily+model.nb_params_cov-1);i++) {
#                res(0,i+1) = -dS1[i];
#                res(i+1,0) = -dS1[i];
#                res(i+1,i+1) = -alpha*d2S1[i*(i+1)/2+i];
#                for(j=0;j<i;j++) {
#                  res(i+1,j+1) = -alpha*d2S1[i*(i+1)/2+j];
#                  res(j+1,i+1) = res(i+1,j+1);
#                }
#             }
#             param[0]=alpha;
#             model.set_params(param);//also memorize the current value for alpha which is not 1 in fact

#         }
#         param[0]=alpha;//LD:changed for bayesian
#         return res;
#     }

#     NumericVector get_alpha_est(NumericVector param) {
#         NumericVector res(1);
#         contrast(param); //To compute S1 and S0
#         res[0]=S0/S1;
#         return res;
#     }

#     VamModel* get_model() {
#     	return model;
#     }


#     //delegate from model cache!
#     List get_virtual_age_infos(double by,double from, double to) {
#         return model.get_virtual_age_infos(by,from,to);
#     }

#     DataFrame get_selected_data(int i) {
#         return model.get_selected_data(i);
#     }

#     void set_covariates(List covariates_) {
#         model.set_covariates(covariates_);
#     }


# function init_mle_vam(mle::MLE;deriv::Bool=false)
#     init!(mle.comp,deriv)
# end

# Rcpp -> init_mle_vam_for_current_system
function init_mle(mle::MLE; deriv::Bool=false)
     
    for mm in mle.model.models
        init!(mm)
    end

    mle.model.Vright = 0
    mle.model.Vright = 0
    mle.model.A = 1
    mle.model.k = 1
    mle.model.id_mod = 0 #id of current model
    init!(mle.model.comp, deriv=deriv)

    for type in mle.model.type
        if type < 0 
            mle.model.comp.S0 += 1
        end
    end
    if deriv
        mle.model.dVright = zeros(mle.model.comp.nbd)
        mle.model.dA = zeros(mle.model.comp.nbd)
        nb2m = mle.model.nb_params_maintenance * (mle.model.nb_params_maintenance + 1) ÷ 2
        mle.model.d2Vright = zeros(nb2m)
        mle.model.d2A = zeros(nb2m)
    end
end

#     void gradient_update_for_current_system() {
#         int i,ii;
#     	contrast_update_for_current_system(true,false);

#         double *cumhVright_param_derivative=model.family.cumulative_hazardRate_param_derivative(model.Vright,true);
#         double *cumhVleft_param_derivative=model.family.cumulative_hazardRate_param_derivative(model.Vleft,false);
#         double *hVleft_param_derivative=model.family.hazardRate_param_derivative(model.Vleft,false);
#         for(i=0;i<model.nb_paramsFamily-1;i++){
#             if(model.k >= left_censor) model.dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i] ;
#             model.dS2[i] += hVleft_param_derivative[i]/model.hVleft*model.indType ;
#         }
#     	double hVright=model.family.hazardRate(model.Vright);
#     	double dhVleft=model.family.hazardRate_derivative(model.Vleft);
#     	  //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model.k,hVright,dhVleft,model.indType);
#     	for(ii=0;ii<model.nb_params_maintenance;ii++,i++) {
#     		if(model.k >= left_censor) model.dS1[i] += model.hVleft * model.dVleft[ii] - hVright * model.dVright[ii];
#     		//printf("dS1[%d]=(%lf,%lf,%lf),%lf,",ii+1,model.hVleft,model.dVleft[ii],model.dVright[ii],model.dS1[ii+1]);
#     		model.dS2[i] +=  dhVleft * model.dVleft[ii]/model.hVleft * model.indType;
#     		//printf("dS2[%d]=%lf,",ii+1,model.dS2[ii+1]);
#             model.dS3[ii] +=  model.dA[ii]/model.A * model.indType;
#     	}
#     	//printf("\n");
#     }

#     void hessian_update_for_current_system() {
#         int i;
#         int j;
#         contrast_update_for_current_system(true,true);

#         double *cumhVright_param_derivative=model.family.cumulative_hazardRate_param_derivative(model.Vright,true);
#         double *cumhVleft_param_derivative=model.family.cumulative_hazardRate_param_derivative(model.Vleft,false);
#         double *hVleft_param_derivative=model.family.hazardRate_param_derivative(model.Vleft,false);
#         double *cumhVright_param_2derivative=model.family.cumulative_hazardRate_param_2derivative(model.Vright,true);
#         double *cumhVleft_param_2derivative=model.family.cumulative_hazardRate_param_2derivative(model.Vleft,false);
#         double *hVleft_param_2derivative=model.family.hazardRate_param_2derivative(model.Vleft);
#         for(i=0;i<model.nb_paramsFamily-1;i++){
#             if(model.k >= left_censor) model.dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i] ;
#             model.dS2[i] += hVleft_param_derivative[i]/model.hVleft*model.indType ;
#             for(j=0;j<=i;j++) {
#                 if(model.k >= left_censor) model.d2S1[i*(i+1)/2+j] += cumhVleft_param_2derivative[i*(i+1)/2+j]-cumhVright_param_2derivative[i*(i+1)/2+j];
#                 model.d2S2[i*(i+1)/2+j] += (hVleft_param_2derivative[i*(i+1)/2+j]/model.hVleft -hVleft_param_derivative[i]*hVleft_param_derivative[j]/pow(model.hVleft,2))*model.indType;
#             }
#         }
#         double hVright=model.family.hazardRate(model.Vright);
#         double dhVleft=model.family.hazardRate_derivative(model.Vleft);
#         double dhVright=model.family.hazardRate_derivative(model.Vright);
#         double *hVright_param_derivative=model.family.hazardRate_param_derivative(model.Vright,true);
#         double *dhVleft_param_derivative=model.family.hazardRate_derivative_param_derivative(model.Vleft);
#         double d2hVleft=model.family.hazardRate_2derivative(model.Vleft);
#         //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model.k,hVright,dhVleft,model.indType);
#         for(i=0;i<model.nb_params_maintenance;i++) {
#             if(model.k >= left_censor) model.dS1[i+model.nb_paramsFamily-1] += model.hVleft * model.dVleft[i] - hVright * model.dVright[i];
#             //printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model.hVleft,model.dVleft[i],model.dVright[i],model.dS1[i+1]);
#             model.dS2[i+model.nb_paramsFamily-1] +=  dhVleft * model.dVleft[i]/model.hVleft * model.indType;
#             //printf("dS2[%d]=%lf,",i+1,model.dS2[i+1]);
#             //column 0 and i+1 corresponds to the line indice of (inferior diagonal part of) the hessian matrice
#             model.dS3[i] +=  model.dA[i]/model.A * model.indType;
#             for(j=0;j<model.nb_paramsFamily-1;j++){
#                 if(model.k >= left_censor) model.d2S1[(i+model.nb_paramsFamily-1)*(i+model.nb_paramsFamily)/2+j] += hVleft_param_derivative[j] * model.dVleft[i] - hVright_param_derivative[j] * model.dVright[i];
#                 model.d2S2[(i+model.nb_paramsFamily-1)*(i+model.nb_paramsFamily)/2+j] +=  dhVleft_param_derivative[j] * model.dVleft[i]/model.hVleft * model.indType - hVleft_param_derivative[j]*dhVleft * model.dVleft[i]/pow(model.hVleft,2) * model.indType;
#             }
#             for(j=0;j<=i;j++){
#                 //i+1 and j+1(<=i+1) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                 if(model.k >= left_censor) model.d2S1[(i+model.nb_paramsFamily-1)*(i+model.nb_paramsFamily)/2+j+model.nb_paramsFamily-1] += dhVleft*model.dVleft[i]*model.dVleft[j] + model.hVleft * model.d2Vleft[i*(i+1)/2+j] - dhVright*model.dVright[i]*model.dVright[j] - hVright * model.d2Vright[i*(i+1)/2+j];
#                 model.d2S2[(i+model.nb_paramsFamily-1)*(i+model.nb_paramsFamily)/2+j+model.nb_paramsFamily-1] += ( model.dVleft[i]*model.dVleft[j]*(d2hVleft/model.hVleft-pow(dhVleft/model.hVleft,2)) + dhVleft * model.d2Vleft[i*(i+1)/2+j]/model.hVleft )* model.indType;
#                 model.d2S3[i*(i+1)/2+j] += (model.d2A[i*(i+1)/2+j]/model.A -model.dA[i]*model.dA[j]/pow(model.A,2))* model.indType;
#             }
#         }
#     }

#     void gradient_dS_maintenance_update(int i,int ii) {
#         dS1[i] += model.dS1[i] * (model.nb_params_cov > 0 ? exp(model.sum_cov) : 1.0); dS2[i] += model.dS2[i]; dS3[ii] += model.dS3[ii];
#     }

#     void gradient_dS_family_update(int i) {
#         dS1[i] += model.dS1[i] * (model.nb_params_cov > 0 ? exp(model.sum_cov) : 1.0); dS2[i] += model.dS2[i];
#     }

#     void gradient_dS_covariate_update(int i,int ii) {
#         //nb_params_cov > 0 necessarily
#         double cov=model.get_covariate(ii);
#         dS1[i] += model.S1 * cov * exp(model.sum_cov);
#         //dS2[i]=0
#         dS4[ii] += model.S0 * cov;
#     }

# };

# #endif //RCPP_MLE_VAM_H
