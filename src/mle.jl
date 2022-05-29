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
    m.left_censor = 0
end

function select_left_censor(mle::MLE, i::Int)
    if length(mle.left_censors) >= i 
        mle.left_censor=left_censors[i]
    end
end

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

function contrast(mle::MLE, param::Vector{Float64}; alpha_fixed::Bool=false)::Float64
    res = 0
    alpha = param[1] #;//save current value of alpha

    param[1] = 1 #//Rmk: alpha replaces param[0] => a bit weird!

    init!(mle.comp)

    params!(mle.model,param);
    # //printf("System %d\n",1);
    select_data(mle.model, 1)
    # if(model.nb_params_cov > 0) mle.model.select_current_system(0,true);
    select_left_censor(mle, 1)
    contrast_current(mle)
    # //only if multi-system
    if mle.model.nb_system > 1
        for i in 2:mle.model.nb_system 
            #//printf("System %d\n",i+1);
            select_data(mle.model, i)
            if mle.model.nb_params_cov > 0 
                select_current_system(mle.model, i, true)
            end
            select_left_censor(mle, i)
            contrast_current(mle)
        end
    end

    # //DEBUG: println("alpha=$alpha,S0=$(mle.comp.S0),S1=$(mle.comp.S1),S2=$(mle.comp.S2),S3=$(mle.comp.S3),S4=$(mle.comp.S4)")
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
    #//model updated for current system: mle.comp.S1,S2,S0,dS1,dS2
    tmp = mle.model.comp.S1
    mle.comp.S2 += mle.model.comp.S2
    mle.comp.S0 += mle.model.comp.S0
    mle.comp.S3 += mle.model.comp.S3
    if mle.model.nb_params_cov > 0
        compute_covariates(mle.model) #//initialize model.sum_cov
        tmp *= exp(mle.model.sum_cov)
        mle.comp.S4 += mle.model.comp.S0 * mle.model.sum_cov
        #//printf("(S0=%lf) * (sum_cov=%lf) = (S4 =%lf)\n",model.S0, mle.model.sum_cov,model.S0 * mle.model.sum_cov);
    end
    mle.comp.S1 += tmp
    #//printf("Conclusion : mle.comp.S1=%f, S2=%f, S0=%f, S4=%f\n",model.S1,model.S2,model.S0,model.S4);
end

function gradient(mle::MLE, param::Vector{Float64}; alpha_fixed::Bool=false)
    res = zeros(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov)
    alpha=param[1] #save current value of alpha

    param[1]=1
    init!(mle.comp, deriv = true)
    params!(mle.model, param)
    select_data(mle.model, 1)
    # if(model.nb_params_cov > 0) mle.model.select_current_system(0,true);
    select_left_censor(mle, 1)
    gradient_current_system(mle)
    # //only if multi-system
    if mle.model.nb_system > 1
        for i in 2:mle.model.nb_system 
            select_data(mle.model, i)
            if mle.model.nb_params_cov > 0 
                select_current_system(mle.model, i, true)
            end
            select_left_censor(mle, i)
            gradient_current(mle)
        end
    end
    # compute gradient
    param[1] = alpha_fixed ? alpha : mle.comp.S0 / mle.Comp.S1

    params!(mle.model, param) # also memorize the current value for alpha which is not 1 in fact

    res[1] = alpha_fixed ? mle.comp.S0/alpha - mle.comp.S1 : 0
    
    np = 1
    for i in 1:mle.model.nb_params_family
        res[i + np] = -mle.comp.dS1[i] * param[1] + mle.comp.dS2[i]
    end
    np += mle.model.nb_params_family
    if mle.model.nb_params_maintenance > 0
        for i in 1:mle.model.nb_params_maintenance
            res[i + np] = -mle.comp.dS1[i + np - 1] * param[1] + mle.comp.dS2[i + np - 1] + mle.comp.dS3[i]
        end
    end
    np += mle.nb_params_maintenance
    if mle.model.nb_params_cov > 0
    for i in 1:model.nb_params_cov
        res[i + np] = -mle.comp.dS1[i + np - 1] * param[1] + mle.comp.dS4[i]
    end

    param[1] = alpha ## BIZARRE!
    return res
end

function gradient_current(mle::MLE)
    init_mle(mle, deriv = true)
    n = length(mle.model.time)
    while mle.model.k < n
        gradient_update_current(mle)
        type = mle.model.type[mle.model.k + 1]
        if type < 0 
            type = 0
        end
        # //model.indMode = (type < 0 ? 0 : type)
        update!(mle.model.models[1 + type], mle.model, gradient = true)
    end
    contrast_update_S(mle)
    #//precomputation of covariate term to multiply (in fact just exp)
    for i = 1:mle.model.nb_params_family
        gradient_update_dS_family(mle, i)
    end
    np = mle.model.nb_params_family
    for i in 1:mle.model.nb_params_maintenance
        gradient_update_dS_maintenance(i + np,i)
    end
    np += mle.model.nb_params_cov
    for i in 1:mle.model.nb_params_cov
        gradient_update_dS_covariate(i + np, i)
    end
end

function gradient_update_current(mle::MLE)
    contrast_update_current(mle, deriv =true)

    cumhVright_param_derivative = cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vright, true)
    cumhVleft_param_derivative=cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    hVleft_param_derivative=hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    for i in 1:mle.model.nb_params_family
        if mle.model.k >= mle.left_censor 
            mle.model.comp.dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i]
        end
        mle.model.comp.dS2[i] += hVleft_param_derivative[i] / mle.model.hVleft * mle.model.indType
    end
    hVright=hazard_rate(mle.model.family, mle.model.Vright)
    dhVleft=hazard_rate_derivative(mle.model.family, mle.model.Vleft)
    # printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model.k,hVright,dhVleft,model.indType);
    np = mle.model.nb_params_family
    for i in 1:mle.model.nb_params_maintenance
        if mle.model.k >= mle.left_censor 
            mle.model.comp.dS1[i + np] += mle.model.hVleft * mle.model.dVleft[i] - hVright * mle.model.dVright[i]
        end
        # printf("dS1[%d]=(%lf,%lf,%lf),%lf,",ii+1,model.hVleft,model.dVleft[ii],model.dVright[ii],model.dS1[ii+1]);
        mle.model.comp.dS2[i + np] +=  dhVleft * mle.model.dVleft[i] / mle.model.hVleft * mle.model.indType
        #//printf("dS2[%d]=%lf,",ii+1,model.dS2[ii+1]);
        mle.model.comp.dS3[i] +=  mle.model.dA[i] / mle.model.A * mle.model.indType
    end
    #//printf("\n");
end

function gradient_update_dS_maintenance(mle::MLE,  i::Int, ii::Int) 
    mle.comp.dS1[i] += mle.model.comp.dS1[i] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
    mle.comp.dS2[i] += mle.model.comp.dS2[i]
    mle.comp.dS3[ii] += mle.model.comp.dS3[ii]
end

function gradient_update_dS_family(mle::MLE, i::Int)
    mle.comp.dS1[i] += mle.model.comp.dS1[i] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
    mle.comp.dS2[i] += mle.model.comp.dS2[i]
end

function gradient_update_dS_covariate(mle::MLE, i::Int, ii::Int)
    #//nb_params_cov > 0 necessarily
    cov=covariate(mle.model, ii)
    mle.comp.dS1[i] += mle.model.comp.S1 * cov * exp(mle.model.sum_cov)
    #//dS2[i]=0
    mle.comp.dS4[ii] += mle.model.comp.S0 * cov
end

function hessian(mle::MLE, param::Vector{Float64}; alpha_fixed::Bool=false)
    res = zeros(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov, mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov)
    alpha=param[1] #save current value of alpha

    param[1]=1
    init!(mle.comp, deriv = true)
    params!(mle.model, param)
    select_data(mle.model, 1)
    # if(model.nb_params_cov > 0) mle.model.select_current_system(0,true);
    select_left_censor(mle, 1)
    hessian_current(mle)

    # //only if multi-system
    if mle.model.nb_system > 1
        for i in 2:mle.model.nb_system 
            select_data(mle.model, i)
            if mle.model.nb_params_cov > 0 
                select_current_system(mle.model, i, true)
            end
            select_left_censor(mle, i)
            hessian_current(mle)
        end
    end

    # //compute hessian
    if !alpha_fixed
        res[1,1] = 0
        for i in 1:(mle.model.nb_params_family - 1)
            res[1, i + 1] = 0
            res[i + 1, 1] = 0
            res[i + 1, i + 1] = mle.comp.dS1[i]^2 / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[i * (i + 1) ÷ 2 + i]/mle.comp.S1 * mle.comp.S0 + mle.comp.d2S2[i * (i + 1) ÷2 + i]
            for j in 1:i # ?? or for j in 0:(i - 1)
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                res[i + 1, j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[i * (i + 1) ÷ 2 + j] / mle.comp.S1 * S0 + d2S2[i * (i + 1) ÷ 2 + j]
                res[j + 1, i + 1] = res[i + 1,j + 1]
            end
        end
        for i in mle.model.nb_params_family:(mle.model.nb_params_maintenance + mle.model.nb_params_family - 1)
            res[1, i + 1] = 0
            res[i + 1, 1] = 0
            res[i + 1, i + 1] = mle.comp.dS1[i]^2 / mle.comp.S1^2 * S0 - d2S1[i * (i + 1) ÷ 2 + i] / mle.comp.S1 * S0 + d2S2[i * (i + 1) ÷ 2 + i] + d2S3[(i-(mle.model.nb_params_family - 1)) * (i -(mle.model.nb_params_family - 1) + 1) ÷ 2 + i - (model.nb_params_family-1)]
            for j in 1:(mle.model.nb_params_family-1)
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                res[i + 1, j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * S0 - d2S1[i * (i + 1) ÷ 2 + j] / mle.comp.S1 * S0 + d2S2[i * (i + 1) ÷ 2 + j]
                res[j + 1, i + 1] = res[i + 1, j + 1]
            end
            for j in mle.model.nb_params_family:(i - 1)
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                res[i + 1, j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * S0 - d2S1[i * (i + 1) ÷ 2 + j] / mle.comp.S1 * S0 + d2S2[i * (i + 1) ÷ 2 + j] + d2S3[(i - (mle.model.nb_params_family-1)) * (i - (mle.model.nb_params_family - 1) + 1) ÷ 2 + j - (mle.model.nb_params_family - 1)]
                res[j + 1, i + 1] = res[i + 1, j + 1]
            end
        end
        for i in (mle.model.nb_params_maintenance + mle.model.nb_params_family):(mle.model.nb_params_maintenance + mle.model.nb_params_family + mle.model.nb_params_cov - 1)
            res[1,i + 1] = 0
            res[i + 1, 1] = 0
            res[i + 1,i + 1] = mle.comp.dS1[i] ^2 / mle.comp.S1^2 * S0 - d2S1[i * (i + 1) ÷ 2 + i] / mle.comp.S1 * S0
            for j in 1:i
                res[i + 1,j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * S0 - d2S1[i * (i + 1) ÷ 2 + j] / mle.comp.S1 * S0
                res[j + 1,i + 1] = res[i + 1,j + 1]
            end
        end
        param[1] = S0/S1
        params!(mle.model, param) #;//also memorize the current value for alpha which is not 1 in fact
    else

        res[1, 1] = -S0 / alpha^2
        for i in 1:(mle.model.nb_params_family-1)
            res[1,i + 1] = -mle.comp.dS1[i]
            res[i + 1,1] = -mle.comp.dS1[i]
            res[i + 1,i + 1] = d2S2[i * (i + 1) ÷ 2 + i] - alpha * d2S1[i *(i + 1) ÷ 2 + i]
            for j in 1:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                res[i + 1,j + 1] = d2S2[i * (i + 1) ÷ 2 + j] - alpha * d2S1[i * (i + 1) ÷2 + j]
                res[j + 1,i + 1] = res[i + 1,j + 1]
            end
        end
        for i in mle.model.nb_params_family:(mle.model.nb_params_maintenance + mle.model.nb_params_family - 1)
            res[1, i + 1] = -mle.comp.dS1[i]
            res[i + 1, 1] = -mle.comp.dS1[i]
            res[i + 1, i + 1] = d2S2[i * (i + 1) ÷ 2 + i] - alpha * d2S1[i * (i + 1) ÷ 2 + i] + d2S3[(i-(mle.model.nb_params_family-1)) * (i - (mle.model.nb_params_family - 1) + 1) ÷ 2 + i - (mle.model.nb_params_family - 1)]
            for j in 1:(mle.model.nb_params_family - 1)
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                res[i + 1, j + 1] = d2S2[i * (i + 1) ÷ 2 + j] - alpha * d2S1[i * (i + 1) ÷ 2 + j]
                res[j + 1, i + 1] = res[i + 1, j + 1]
            end
            for j in mle.model.nb_params_family:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                res[i + 1, j + 1] = d2S2[i * (i + 1) ÷ 2 + j] - alpha * d2S1[i * (i + 1) ÷ 2 + j] + d2S3[(i - (mle.model.nb_params_family - 1)) * (i - (mle.model.nb_params_family - 1) + 1) ÷ 2 + j - (mle.model.nb_params_family - 1)]
                res[j + 1, i + 1] = res[i + 1, j + 1]
            end
        end
        for i in (model.nb_params_maintenance+mle.model.nb_params_family):(mle.model.nb_params_maintenance+mle.model.nb_params_family + mle.model.nb_params_cov - 1)
            res[1, i + 1] = -mle.comp.dS1[i]
            res[i + 1, 1] = -mle.comp.dS1[i]
            res[i + 1, i + 1] = -alpha * d2S1[i * (i + 1) ÷ 2 + i]
            for j in 1:i
                res[i + 1, j + 1] = -alpha * d2S1[i * (i + 1) ÷ 2 + j]
                res[j + 1, i + 1] = res[i + 1 , j + 1]
            end
        end
        param[1] = alpha
        params!(mle.model, param) #also memorize the current value for alpha which is not 1 in fact

    end
    param[1] = alpha # LD:changed for bayesian
    return res
end

function hessian_current(mle::MLE)
    init_mle(mle, deriv = true)
    n = length(mle.model.time)
    while mle.model.k < n
        hessian_update_current(mle)
        type = mle.model.type[mle.model.k + 1]
        if type < 0 
            type = 0
        end
        # //model.indMode = (type < 0 ? 0 : type)
        update!(mle.model.models[1 + type], mle.model, hessian = true)
    end
    contrast_update_S(mle)
    for i = 1:mle.model.nb_params_family
        gradient_update_dS_family(mle, i)
        for j in 1:(i + 1)
            #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            k = i * (i + 1) ÷ 2 + j
            d2S1[k] += mle.model.d2S1[k] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
            d2S2[k] += mle.model.d2S2[k]
        end
    end
    np = mle.model.nb_params_family
    for i in 1:mle.model.nb_params_maintenance
        gradient_update_dS_maintenance(mle, i + np,i)
        for j in 1:(i + 1)
            #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            k = (i + np) * (i + np + 1) ÷ 2 + j
            d2S1[k] += mle.model.d2S1[k] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
            d2S2[k] += mle.model.d2S2[k]
            #//ii and j(<=ii) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            kk = i * (i + 1) ÷ 2 + j
            d2S3[kk] += mle.model.d2S3[kk]
        end
        for j=(i+2):(i + np)
            # //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            k = (i + np) * (i + np +1) ÷ 2 + j
            d2S1[k] += mle.model.d2S1[k] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
            d2S2[k] += mle.model.d2S2[k]
        end
    end
    np += mle.model.nb_params_cov
    for i in 1:mle.model.nb_params_cov
        gradient_update_dS_covariate(mle, i + np, i)
        for j in 1:(model.nb_params_family-1 + mle.model.nb_params_maintenance)
            k = (i + np) * (i + np + 1) ÷ 2 + j
            d2S1[k] += covariate(mle.model, i) * exp(mle.model.sum_cov) * mle.model.dS1[j]
        end
        for j in (mle.model.nb_params_family + mle.model.nb_params_maintenance):(i+1)
            k = (i + np) * (i + np + 1) ÷ 2 + j
            d2S1[k] += covariate(mle.model, i) * covariate(mle.model, j - mle.model.nb_params_family + 1 - mle.model.nb_params_maintenance) * exp(mle.model.sum_cov) * mle.model.S1;
        end
    end
end

function hessian_update_current(mle::MLE) 
    contrast_update_current(mle, deriv =true)

    cumhVright_param_derivative = cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vright, true)
    cumhVleft_param_derivative = cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    hVleft_param_derivative = hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    cumhVright_param_2derivative = cumulative_hazard_rate_param_2derivative(mle.model.family, mle.model.Vright,true)
    cumhVleft_param_2derivative = cumulative_hazard_rate_param_2derivative(mle.model.family, mle.model.Vleft,false)
    hVleft_param_2derivative = hazard_rate_param_2derivative(mle.model.family, mle.model.Vleft)
    for i in 1:(model.nb_params_family - 1)
        if mle.model.k >= mle.left_censor 
            mle.model.comp.dS1[i] +=  cumhVleft_param_derivative[i] - cumhVright_param_derivative[i]
        end
        mle.model.comp.dS2[i] += hVleft_param_derivative[i] / mle.model.hVleft * mle.model.indType ;
        for j in 1:(i + 1)
            ij = i * (i + 1) ÷ 2 + j
            if mle.model.k >= mle.left_censor 
                mle.model.comp.d2S1[ij] += cumhVleft_param_2derivative[ij] - cumhVright_param_2derivative[ij]
            end
            mle.model.comp.d2S2[ij] += (hVleft_param_2derivative[ij] / mle.model.hVleft - hVleft_param_derivative[i] * hVleft_param_derivative[j] / mle.model.hVleft^2) * mle.model.indType
        end
    end
    hVright = hazard_rate(mle.model.family, mle.model.Vright)
    dhVleft = hazard_rate_derivative(mle.model.family, mle.model.Vleft)
    dhVright = hazard_rate_derivative(mle.model.family, mle.model.Vright)
    hVright_param_derivative = hazard_rate_param_derivative(mle.model.family, mle.model.Vright,true)
    dhVleft_param_derivative = hazard_rate_derivative_param_derivative(mle.model.family, mle.model.Vleft)
    d2hVleft = hazard_rate_2derivative(mle.model.family, mle.model.Vleft)
    # //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model.k,hVright,dhVleft,model.indType);
    npf = mle.model.nb_params_family - 1
    npf = mle.model.nb_params_maintenance
    for i in 1:npm
        if mle.model.k >= mle.left_censor
            mle.model.comp.dS1[i + npf] += mle.model.hVleft * mle.model.dVleft[i] - hVright * mle.model.dVright[i];
        end
        # //printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model.hVleft,model.dVleft[i],model.dVright[i],model.dS1[i+1]);
        mle.model.comp.dS2[i + npf] +=  dhVleft * mle.model.dVleft[i]/model.hVleft * mle.model.indType
        #//printf("dS2[%d]=%lf,",i+1,model.dS2[i+1]);
        #//column 0 and i+1 corresponds to the line indice of (inferior diagonal part of) the hessian matrice
        mle.model.comp.dS3[i] +=  mle.model.dA[i]/model.A * mle.model.indType;
        for j in 1:npf
            ij = (i + npf) * (i + npf + 1) ÷ 2 + j
            if mle.model.k >= mle.left_censor
                mle.model.comp.d2S1[ij] += hVleft_param_derivative[j] * mle.model.dVleft[i] - hVright_param_derivative[j] * mle.model.dVright[i]
            end
            mle.model.comp.d2S2[ij] +=  dhVleft_param_derivative[j] * mle.model.dVleft[i] / mle.model.hVleft * mle.model.indType - hVleft_param_derivative[j] * dhVleft * mle.model.dVleft[i] / mle.model.hVleft^2 * mle.model.indType
        end
        for j=1:(i + 1)
            #//i+1 and j+1(<=i+1) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ij = (i + npf) * (i + npf + 1) ÷ 2 + j + npf
            iij = i * (i + 1) ÷ 2 + j
            if mle.model.k >= mle.left_censor
                mle.model.comp.d2S1[ij] += dhVleft*model.dVleft[i] * mle.model.dVleft[j] + mle.model.hVleft * mle.model.d2Vleft[iij] - dhVright * mle.model.dVright[i] * mle.model.dVright[j] - hVright * mle.model.d2Vright[iij]
            end
            mle.model.comp.d2S2[ij] += (mle.model.dVleft[i] * mle.model.dVleft[j] * (d2hVleft / mle.model.hVleft - (dhVleft/model.hVleft)^2) + dhVleft * mle.model.d2Vleft[iij] / mle.model.hVleft) * mle.model.indType
            mle.model.comp.d2S3[iij] += (mle.model.d2A[iij] / mle.model.A - mle.model.dA[i] * mle.model.dA[j] / mle.model.A^2) * mle.model.indType
        end
    end
end

function alpha_est(mle::MLE, param::Vector{Float64})
    contrast(mle, param) #//To compute mle.comp.S1 and S0
    return mle.comp.S0 / mle.comp.S1
end



#     //delegate from model cache!
#     List get_virtual_age_infos(double by,double from, double to) {
#         return mle.model.get_virtual_age_infos(by,from,to);
#     }

#     DataFrame get_selected_data(int i) {
#         return mle.model.get_selected_data(i);
#     }

#     void set_covariates(List covariates_) {
#         mle.model.set_covariates(covariates_);
#     }

