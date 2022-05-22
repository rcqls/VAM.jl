mutable struct MaintenanceModel
    id::Int
    id_params::Int
end

abstract type AbstractMaintenanceModel end
mutable struct ARA1 <: AbstractMaintenanceModel
    ρ::Float64
end
mutable struct ARAInf <: AbstractMaintenanceModel 
	ρ::Float64
end

struct AGAN <: AbstractMaintenanceModel
end
struct ABAO <: AbstractMaintenanceModel
end

struct AGAP <: AbstractMaintenanceModel
end
mutable struct QAGAN <: AbstractMaintenanceModel
    ρ::Float64
end

mutable struct QR <: AbstractMaintenanceModel
    ρ::Float64
end

abstract type F_GQR end
struct ID_GQR <:  F_GQR
end

struct LOG_GQR <: F_GQR
end

struct SQRT_GQR <: F_GQR
end

mutable struct GQR <: AbstractMaintenanceModel
    ρ::Float64
    K::Float64
    f::F_GQR
end
mutable struct GQR_ARA1 <:  AbstractMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
end
mutable struct GQR_ARAInf <:  AbstractMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
end
mutable struct ARAm <: AbstractMaintenanceModel
    ρ::Float64
    m::Int
end
mutable struct GQR_ARAm <: AbstractMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
    m::Int
end

function update(m::ARA1, model::AbstractModel; with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # double prov;
    model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu;
    end
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    if with_hessian
        for k in nk-1:-1:1 
            for i in 0:(model.nb_paramsMaintenance - 1)
                for j in 0:i
                    model.d2VR_prec[k * (model.nb_paramsMaintenance * (model.nb_paramsMaintenance + 1) / 2) + i * (i + 1) / 2 + j ] = model.d2VR_prec[ (k - 1 ) * (model.nb_paramsMaintenance * (model.nb_paramsMaintenance + 1) / 2) + i * (i + 1) / 2 + j]
                end
            end
        end
        if nk > 0 
            for i in 0:(model.nb_paramsMaintenance-1)
                for j in 0:i
                    prov = (1-m.ρ)*model.d2A[i * (i + 1) / 2 + j] * (model.time[model.k]-model.time[model.k - 1])
                    model.d2VR_prec[i * (i + 1) / 2 + j] = prov
                    model.d2Vright[i * (i + 1) / 2 + j] += prov
                end
            end
            for j in 0:id_params
                prov = model.dA[j] * (model.time[model.k] - model.time[model.k - 1])
                model.d2VR_prec[id_params * (id_params + 1) / 2 + j] -= prov
                model.d2Vright[id_params * (id_params + 1) / 2 + j] -= prov
            end
            for i in id_params:(model.nb_paramsMaintenance - 1)
                prov = model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
                model.d2VR_prec[i * (i + 1) / 2 + id_params] -= prov
                model.d2Vright[i * (i + 1) / 2 + id_params] -= prov
            end
        else
           for i in 0:(model.nb_paramsMaintenance-1)
                for j in 0:i
                    model.d2Vright[i * (i + 1) / 2 + j] += (1-m.ρ) * model.d2A[i * (i + 1) / 2 + j] * (model.time[model.k]-model.time[model.k - 1])
                end
            end
            for j in 0:id_params
                model.d2Vright[id_params * (id_params + 1) / 2 + j] -= model.dA[j] * (model.time[model.k] - model.time[model.k - 1])
            end
            for i in id_params:(model.nb_paramsMaintenance - 1)
                model.d2Vright[i * (i + 1) / 2 + id_params] -= model.dA[i]*(model.time[model.k]-model.time[model.k - 1])
            end 
        end
    end
    if with_gradient || with_hessian
        for k in (nk-1):-1:1
            for i in 0:(model.nb_paramsMaintenance - 1)
                model.dVR_prec[k * model.nb_paramsMaintenance + i] = model.dVR_prec[(k-1) * model.nb_paramsMaintenance + i]
            end
        end
        if nk > 0
            for i in 0:(model.nb_paramsMaintenance - 1)
                prov = (1 - m.ρ) * model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
            prov=model.A * (model.time[model.k] - model.time[model.k - 1])
            model.dVR_prec[id_params] -= prov
            model.dVright[id_params] -= prov
        else
            for i in 0:(model.nb_paramsMaintenance - 1)
                model.dVright[i] += (1-m.ρ) * model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
            end
            model.dVright[id_params] -= model.A * (model.time[model.k] - model.time[model.k - 1])
        end
    end
    for k in (nk-1):-1:1
        model.VR_prec[k] = model.VR_prec[k-1]
    end
    prov = (1 - m.ρ) * model.A * (model.time[model.k] - model.time[model.k - 1])
    if nk>0
        model.VR_prec[0] = prov
    end
    model.Vright += prov

    #// save old model
    model.idMod = id
end

function update(m::ARAInf, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # model->k += 1;


    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    # if (with_hessian){
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params;j++) {
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
    #         }
    #         for(i=id_params;i<model->nb_paramsMaintenance;i++) {
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
    #         }
    #     }
    #     if(nk>0){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[i*(i+1)/2+j] = (1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #                 model->d2Vright[i*(i+1)/2+j] = (1-rho)*model->d2Vleft[i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params;j++) {
    #             model->d2VR_prec[id_params*(id_params+1)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2Vright[id_params*(id_params+1)/2+j] -= model->dVleft[j];
    #         }
    #         for(i=id_params;i<model->nb_paramsMaintenance;i++) {
    #             model->d2VR_prec[i*(i+1)/2+id_params] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2Vright[i*(i+1)/2+id_params] -= model->dVleft[i];
    #         }
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j] = (1-rho)*model->d2Vleft[i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params;j++) {
    #             model->d2Vright[id_params*(id_params+1)/2+j] -= model->dVleft[j];
    #         }
    #         for(i=id_params;i<model->nb_paramsMaintenance;i++) {
    #             model->d2Vright[i*(i+1)/2+id_params] -= model->dVleft[i];
    #         }
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #         model->dVR_prec[k*model->nb_paramsMaintenance+id_params]-=model->VR_prec[k-1];
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[i]=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->dVright[i]=(1-rho)*model->dVleft[i];
    #         }
    #         model->dVR_prec[id_params]-= model->A*(model->time[model->k]-model->time[model->k - 1]);
    #         model->dVright[id_params]-= model->Vleft;
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]=(1-rho)*model->dVleft[i];
    #         }
    #         model->dVright[id_params]-= model->Vleft;
    #     }
    # }
    # model->Vright=(1-rho)*model->Vleft;
    # for(k=nk-1;k>0;k--) model->VR_prec[k]=(1-rho)*model->VR_prec[k-1];
    # if(nk>0) model->VR_prec[0]=(1-rho)*model->A*(model->time[model->k]-model->time[model->k - 1]);  
    # // save old model
    # model->idMod = id;
end


function update(m::AGAN, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # model->k += 1;

    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # model->A=1;
    # model->Vright=0;
    # for(k=0;k<nk;k++) {
    #     model->VR_prec[k]=0;
    # }
    # if (with_hessian){
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dA[i] = 0;
    #         model->dVright[i] = 0;
    #         for(k=0;k<nk;k++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
    #         }
    #         for(j=0;j<=i;j++) {
    #             //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #             model->d2A[i*(i+1)/2+j] = 0;
    #             model->d2Vright[i*(i+1)/2+j] = 0;
    #             for(k=0;k<nk;k++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
    #             }
    #         }
    #     }
    # }
    # if(with_gradient) {
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dA[i] = 0;
    #         model->dVright[i] = 0;
    #         for(k=0;k<nk;k++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
    #         }
    #     }
    # }

    # // save old model
    # model->idMod = id;

    # //init QR and GQR type models
    # for(i=0;i<model->nbPM + 1;i++) model->models->at(i)->init();
end

function update(m::ABAO, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # double prov;
    # model->k += 1;

    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # if (with_hessian){
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 prov= model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #                 model->d2VR_prec[i*(i+1)/2+j] = prov;
    #                 model->d2Vright[i*(i+1)/2+j]+=prov;
    #             }
    #         }
    #     } else {
    #        for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j]+=model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #             }
    #         }
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->dVR_prec[i]=prov;
    #             model->dVright[i]+=prov;
    #         }
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]+=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #         }
    #     }
    # }
    # for(k=nk-1;k>0;k--) model->VR_prec[k]=model->VR_prec[k-1];
    # prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
    # if(nk>0) model->VR_prec[0]=prov;
    # model->Vright+=prov;    

    # // save old model
    # model->idMod = id;
end

function update(m::AGAP, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # model->k += 1;

    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # if (with_hessian){
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) model->d2VR_prec[i*(i+1)/2+j] = 0;
    #         }
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) model->dVR_prec[i]=0;
    #     }
    # }
    # for(k=nk-1;k>0;k--) model->VR_prec[k]=model->VR_prec[k-1];
    # if(nk>0) model->VR_prec[0]=0;

    # // save old model
    # model->idMod = id;
end


function update(m::QAGAN, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # model->k += 1;

    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # model->Vright=0;
    # for(k=0;k<nk;k++) {
    #     model->VR_prec[k]=0;
    # }
    # if (with_hessian){
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dVright[i] = 0;
    #         for(k=0;k<nk;k++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
    #         }
    #         for(j=0;j<=i;j++) {
    #             //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #             model->d2Vright[i*(i+1)/2+j] = 0;
    #             for(k=0;k<nk;k++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
    #             }
    #         }
    #     }
    # }
    # if(with_gradient) {
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dVright[i] = 0;
    #         for(k=0;k<nk;k++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
    #         }
    #     }
    # }

    # // save old model
    # model->idMod = id;
end


function update(m::QR, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
# void QR::update(bool with_gradient,bool with_hessian) {
#     int i;
#     int j;
#     int k;
#     model->k += 1;

#     int nk=model->k;
#     if(nk>model->mu){
#         nk=model->mu;
#     }
#     if (with_hessian){
#         for(i=0;i<model->nb_paramsMaintenance;i++) {
#             for(j=0;j<=i;j++) {
#                 //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                 model->d2A[i*(i+1)/2+j] = rho* model->d2A[i*(i+1)/2+j];
#                 model->d2Vright[i*(i+1)/2+j] = 0;
#                 for(k=0;k<nk;k++) {
#                     model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
#                 }
#             }
#         }
#         for(j=0;j<=id_params;j++) {
#             //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
#             model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + model->dA[j];
#         }
#         for(i=id_params;i<model->nb_paramsMaintenance;i++) {
#              //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#             model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + model->dA[i];
#         }
#     }
#     if(with_gradient||with_hessian) {
#         for(i=0;i<model->nb_paramsMaintenance;i++) {
#             model->dA[i] = rho *  model->dA[i];
#             model->dVright[i] = 0;
#             for(k=0;k<nk;k++) {
#                 model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
#             }
#         }
#         model->dA[id_params] = model->dA[id_params] +  model->A;
#     }
#     model->A=rho*model->A;
#     model->Vright=0;
#     for(k=0;k<nk;k++) {
#         model->VR_prec[k]=0;
#     }
#     model->idMod = id;
end

function update(m::GQR, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # model->k += 1;
    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # K++;
    # if (with_hessian){
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         for(j=0;j<=i;j++) {
    #             //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #             model->d2A[i*(i+1)/2+j] = pow(rho,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
    #             model->d2Vright[i*(i+1)/2+j] = 0;
    #             for(k=0;k<nk;k++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=0;
    #             }
    #         }
    #     }
    #     for(j=0;j<id_params;j++) {
    #         //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho*model->dA[j];
    #     }
    #     model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho*model->A);
    #     for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #          //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho*model->dA[i];
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dA[i] = pow(rho,f->eval(K)-f->eval(K-1)) *  model->dA[i];
    #         model->dVright[i] = 0;
    #         for(k=0;k<nk;k++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=0;
    #         }
    #     }
    #     model->dA[id_params] = model->dA[id_params] + pow(rho,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho * model->A;
    # }
    # model->A=pow(rho,f->eval(K)-f->eval(K-1))*model->A;
    # model->Vright=0;
    # for(k=0;k<nk;k++) {
    #     model->VR_prec[k]=0;
    # }
    # model->idMod = id;
end

function update(m::GQR_ARA1, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # double prov;
    # K++;
    # model->k += 1;

    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    # if (with_hessian){
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 prov= (1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #                 model->d2VR_prec[i*(i+1)/2+j] = prov;
    #                 model->d2Vright[i*(i+1)/2+j]+=prov;
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             prov=model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2VR_prec[(id_params+1)*(id_params+2)/2+j] -= prov;
    #             model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= prov;
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2VR_prec[i*(i+1)/2+id_params+1] -= prov;
    #             model->d2Vright[i*(i+1)/2+id_params+1] -= prov;
    #         }
    #     } else {
    #        for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j]+=(1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             model->d2Vright[i*(i+1)/2+id_params+1] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #         } 
    #     }
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         for(j=0;j<=i;j++) {
    #             //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #             model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
    #         }
    #     }
    #     for(j=0;j<id_params;j++) {
    #         //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
    #     }
    #     model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
    #     for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #          //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             prov=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->dVR_prec[i]=prov;
    #             model->dVright[i]+=prov;
    #         }
    #         prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
    #         model->dVR_prec[id_params+1]-= prov;
    #         model->dVright[id_params+1]-=prov;
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]+=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #         }
    #         model->dVright[id_params+1]-=model->A*(model->time[model->k]-model->time[model->k - 1]);
    #     }
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
    #     }
    #     model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;

    # }
    # for(k=nk-1;k>0;k--) model->VR_prec[k]=model->VR_prec[k-1];
    # prov=(1-rho_ARA)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    # if(nk>0) model->VR_prec[0]=prov;
    # model->Vright+=prov;
    # model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;

    # // save old model
    # model->idMod = id;
end

function update(m::GQR_ARAInf, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # K++;
    # model->k += 1;


    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    # if (with_hessian){
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
    #         }
    #     }
    #     if(nk>0){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[i*(i+1)/2+j] = (1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #                 model->d2Vright[i*(i+1)/2+j] = (1-rho_ARA)*model->d2Vleft[i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             model->d2VR_prec[(id_params+1)*(id_params+2)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dVleft[j];
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             model->d2VR_prec[i*(i+1)/2+id_params+1] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVleft[i];
    #         }
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j] = (1-rho_ARA)*model->d2Vleft[i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dVleft[j];
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVleft[i];
    #         }
    #     }
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         for(j=0;j<=i;j++) {
    #             //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #             model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
    #         }
    #     }
    #     for(j=0;j<id_params;j++) {
    #         //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
    #     }
    #     model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
    #     for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #          //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(k=nk-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #         model->dVR_prec[k*model->nb_paramsMaintenance+id_params+1]-=model->VR_prec[k-1];
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[i]=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->dVright[i]=(1-rho_ARA)*model->dVleft[i];
    #         }
    #         model->dVR_prec[id_params+1]-= model->A*(model->time[model->k]-model->time[model->k - 1]);
    #         model->dVright[id_params+1]-= model->Vleft;
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]=(1-rho_ARA)*model->dVleft[i];
    #         }
    #         model->dVright[id_params+1]-= model->Vleft;
    #     }
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
    #     }
    #     model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
    # }
    # model->Vright=(1-rho_ARA)*model->Vleft;
    # for(k=nk-1;k>0;k--) model->VR_prec[k]=(1-rho_ARA)*model->VR_prec[k-1];
    # if(nk>0) model->VR_prec[0]=(1-rho_ARA)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    # model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;

    # // save old model
    # model->idMod = id;
end

function update(m::ARAm, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # double prov;
    # model->k += 1;


    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # int nk2=nk;
    # if(nk>m-1){
    #     nk2=m-1;
    # }

    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    # if (with_hessian){
    #     for(k=nk-1;k>nk2;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #     }

    #     if ((model->k>=m)&&(nk2>0)) {
    #         if(nk>nk2) {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) {
    #                 for(j=0;j<=i;j++) {
    #                     model->d2Vright[i*(i+1)/2+j]-=rho*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #                     model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #                 }
    #             }
    #             for(j=0;j<=id_params;j++) {
    #                 model->d2Vright[id_params*(id_params+1)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
    #                 model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
    #             }
    #             for(i=id_params;i<model->nb_paramsMaintenance;i++) {
    #                 model->d2Vright[i*(i+1)/2+id_params] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
    #                 model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
    #             }
    #         } else {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) {
    #                 for(j=0;j<=i;j++) model->d2Vright[i*(i+1)/2+j]-=rho*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #             for(j=0;j<=id_params;j++) model->d2Vright[id_params*(id_params+1)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
    #             for(i=id_params;i<model->nb_paramsMaintenance;i++) model->d2Vright[i*(i+1)/2+id_params] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
    #         }
    #     }

    #     for(k=nk2-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j]-=rho*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params;j++) {
    #             model->d2Vright[id_params*(id_params+1)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+id_params*(id_params+1)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
    #         }
    #         for(i=id_params;i<model->nb_paramsMaintenance;i++) {
    #             model->d2Vright[i*(i+1)/2+id_params] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 prov= (1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #                 model->d2VR_prec[i*(i+1)/2+j] = prov;
    #                 model->d2Vright[i*(i+1)/2+j]+=prov;
    #             }
    #         }
    #         for(j=0;j<=id_params;j++) {
    #             prov=model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2VR_prec[id_params*(id_params+1)/2+j] -= prov;
    #             model->d2Vright[id_params*(id_params+1)/2+j] -= prov;
    #         }
    #         for(i=id_params;i<model->nb_paramsMaintenance;i++) {
    #             prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2VR_prec[i*(i+1)/2+id_params] -= prov;
    #             model->d2Vright[i*(i+1)/2+id_params] -= prov;
    #         }
    #     } else {
    #        for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j]+=(1-rho)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #             }
    #         }
    #         for(j=0;j<=id_params;j++) {
    #             model->d2Vright[id_params*(id_params+1)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #         }
    #         for(i=id_params;i<model->nb_paramsMaintenance;i++) {
    #             model->d2Vright[i*(i+1)/2+id_params] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #         } 
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(k=nk-1;k>nk2;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #     }
    #     if ((model->k>=m)&&(nk2>0)) {
    #         if(nk>nk2) {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) {
    #                 model->dVright[i]-=rho*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
    #                 model->dVR_prec[nk2*model->nb_paramsMaintenance+i]=(1-rho)*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
    #             }
    #             model->dVR_prec[nk2 * (model->nb_paramsMaintenance)+id_params] -= model->VR_prec[nk2-1];
    #             model->dVright[id_params]-=model->VR_prec[nk2-1];
    #         } else {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) model->dVright[i]-=rho*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
    #             model->dVright[id_params]-=model->VR_prec[nk2-1];
    #         }
    #     }
    #     for(k=nk2-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]-=rho*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #         model->dVR_prec[k * (model->nb_paramsMaintenance)+id_params] -= model->VR_prec[k-1];
    #         model->dVright[id_params]-=model->VR_prec[k-1];
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             prov=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->dVR_prec[i]=prov;
    #             model->dVright[i]+=prov;
    #         }
    #         prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
    #         model->dVR_prec[id_params]-= prov;
    #         model->dVright[id_params]-=prov;
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]+=(1-rho)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #         }
    #         model->dVright[id_params]-=model->A*(model->time[model->k]-model->time[model->k - 1]);
    #     }
    # }
    # //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 
    # for(k=nk-1;k>nk2;k--) {
    #     model->VR_prec[k]=model->VR_prec[k-1];
    # }
    # if ((model->k>=m)&&(nk2>0)) {
    #     if(nk>nk2) {
    #         model->Vright-=rho*model->VR_prec[nk2-1];
    #         model->VR_prec[nk2]=(1-rho)*model->VR_prec[nk2-1];
    #     } else model->Vright-=rho*model->VR_prec[nk2-1];
    # }
    # for(k=nk2-1;k>0;k--) {
    #     model->Vright-=rho*model->VR_prec[k-1];
    #     model->VR_prec[k]=(1-rho)*model->VR_prec[k-1];
    # } 
    # prov=(1-rho)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    # //printf("Vright=%f, rho=%f, A=%f, Tk=%f, Tk-1=%f\n",model->Vright,rho,model->A,model->time[model->k],model->time[model->k - 1]);
    # if(nk>0) model->VR_prec[0]=prov;
    # model->Vright+=prov;
    # //printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model->Vright,nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 

    # // save old model
    # model->idMod = id;
end

function update(m::GQR_ARAm, model::AbstractModel;with_gradient::Bool=false,with_hessian::Bool=false)
    # int i;
    # int j;
    # int k;
    # double prov;
    # model->k += 1;
    # K++;

    # int nk=model->k;
    # if(nk>model->mu){
    #     nk=model->mu;
    # }
    # int nk2=nk;
    # if(nk>m-1){
    #     nk2=m-1;
    # }
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model->k,model->max_mem,nk);
    # if (with_hessian){
    #     for(k=nk-1;k>nk2;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #     }

    #     if ((model->k>=m)&&(nk2>0)) {
    #         if(nk>nk2) {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) {
    #                 for(j=0;j<=i;j++) {
    #                     model->d2Vright[i*(i+1)/2+j]-=rho_ARA*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #                     model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #                 }
    #             }
    #             for(j=0;j<=id_params+1;j++) {
    #                 model->d2Vright[(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
    #                 model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
    #             }
    #             for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #                 model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
    #                 model->d2VR_prec[nk2*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
    #             }
    #         } else {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) {
    #                 for(j=0;j<=i;j++) model->d2Vright[i*(i+1)/2+j]-=rho_ARA*model->d2VR_prec[(nk2-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #             for(j=0;j<=id_params+1;j++) model->d2Vright[(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+j];
    #             for(i=id_params+1;i<model->nb_paramsMaintenance;i++) model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVR_prec[(nk2-1)*(model->nb_paramsMaintenance)+i];
    #         }
    #     }

    #     for(k=nk2-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j]-=rho_ARA*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #                 model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j]=(1-rho_ARA)*model->d2VR_prec[(k-1)*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+j];
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             model->d2Vright[(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+(id_params+1)*(id_params+2)/2+j]-= model->dVR_prec[(k-1)*model->nb_paramsMaintenance+j];
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             model->d2Vright[i*(i+1)/2+id_params+1] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
    #             model->d2VR_prec[k*(model->nb_paramsMaintenance*(model->nb_paramsMaintenance+1)/2)+i*(i+1)/2+id_params+1] -= model->dVR_prec[(k-1)*(model->nb_paramsMaintenance)+i];
    #         }
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 prov= (1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #                 model->d2VR_prec[i*(i+1)/2+j] = prov;
    #                 model->d2Vright[i*(i+1)/2+j]+=prov;
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             prov=model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2VR_prec[(id_params+1)*(id_params+2)/2+j] -= prov;
    #             model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= prov;
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             prov=model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->d2VR_prec[i*(i+1)/2+id_params+1] -= prov;
    #             model->d2Vright[i*(i+1)/2+id_params+1] -= prov;
    #         }
    #     } else {
    #        for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             for(j=0;j<=i;j++) {
    #                 model->d2Vright[i*(i+1)/2+j]+=(1-rho_ARA)*model->d2A[i*(i+1)/2+j]*(model->time[model->k]-model->time[model->k - 1]);
    #             }
    #         }
    #         for(j=0;j<=id_params+1;j++) {
    #             model->d2Vright[(id_params+1)*(id_params+2)/2+j] -= model->dA[j]*(model->time[model->k]-model->time[model->k - 1]);
    #         }
    #         for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #             model->d2Vright[i*(i+1)/2+id_params+1] -= model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #         } 
    #     }
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         for(j=0;j<=i;j++) {
    #             //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #             model->d2A[i*(i+1)/2+j] = pow(rho_QR,f->eval(K)-f->eval(K-1))* model->d2A[i*(i+1)/2+j];
    #         }
    #     }
    #     for(j=0;j<id_params;j++) {
    #         //i(<=id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[id_params*(id_params+1)/2+j] = model->d2A[id_params*(id_params+1)/2+j] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[j];
    #     }
    #     model->d2A[id_params*(id_params+1)/2+id_params] = model->d2A[id_params*(id_params+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*(2*model->dA[id_params]+(f->eval(K)-f->eval(K-1)-1)/rho_QR*model->A);
    #     for(i=id_params+1;i<model->nb_paramsMaintenance;i++) {
    #          //id and i(>=id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
    #         model->d2A[i*(i+1)/2+id_params] = model->d2A[i*(i+1)/2+id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1))*(f->eval(K)-f->eval(K-1))/rho_QR*model->dA[i];
    #     }
    # }
    # if(with_gradient||with_hessian) {
    #     for(k=nk-1;k>nk2;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #     }
    #     if ((model->k>=m)&&(nk2>0)) {
    #         if(nk>nk2) {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) {
    #                 model->dVright[i]-=rho_ARA*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
    #                 model->dVR_prec[nk2*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
    #             }
    #             model->dVR_prec[nk2 * (model->nb_paramsMaintenance)+id_params+1] -= model->VR_prec[nk2-1];
    #             model->dVright[id_params+1]-=model->VR_prec[nk2-1];
    #         } else {
    #             for(i=0;i<model->nb_paramsMaintenance;i++) model->dVright[i]-=rho_ARA*model->dVR_prec[(nk2-1)*model->nb_paramsMaintenance+i];
    #             model->dVright[id_params+1]-=model->VR_prec[nk2-1];
    #         }
    #     }
    #     for(k=nk2-1;k>0;k--){
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]-=rho_ARA*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #             model->dVR_prec[k*model->nb_paramsMaintenance+i]=(1-rho_ARA)*model->dVR_prec[(k-1)*model->nb_paramsMaintenance+i];
    #         }
    #         model->dVR_prec[k * (model->nb_paramsMaintenance)+id_params+1] -= model->VR_prec[k-1];
    #         model->dVright[id_params+1]-=model->VR_prec[k-1];
    #     }
    #     if(nk>0) {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             prov=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #             model->dVR_prec[i]=prov;
    #             model->dVright[i]+=prov;
    #         }
    #         prov=model->A*(model->time[model->k]-model->time[model->k - 1]);
    #         model->dVR_prec[id_params+1]-= prov;
    #         model->dVright[id_params+1]-=prov;
    #     } else {
    #         for(i=0;i<model->nb_paramsMaintenance;i++) {
    #             model->dVright[i]+=(1-rho_ARA)*model->dA[i]*(model->time[model->k]-model->time[model->k - 1]);
    #         }
    #         model->dVright[id_params+1]-=model->A*(model->time[model->k]-model->time[model->k - 1]);
    #     }
    #     for(i=0;i<model->nb_paramsMaintenance;i++) {
    #         model->dA[i] = pow(rho_QR,f->eval(K)-f->eval(K-1)) *  model->dA[i];
    #     }
    #     model->dA[id_params] = model->dA[id_params] + pow(rho_QR,f->eval(K)-f->eval(K-1)) * (f->eval(K)-f->eval(K-1))/rho_QR * model->A;
    # }
    # //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 
    # for(k=nk-1;k>nk2;k--) {
    #     model->VR_prec[k]=model->VR_prec[k-1];
    # }
    # if ((model->k>=m)&&(nk2>0)) {
    #     if(nk>nk2) {
    #         model->Vright-=rho_ARA*model->VR_prec[nk2-1];
    #         model->VR_prec[nk2]=(1-rho_ARA)*model->VR_prec[nk2-1];
    #     } else model->Vright-=rho_ARA*model->VR_prec[nk2-1];
    # }
    # for(k=nk2-1;k>0;k--) {
    #     model->Vright-=rho_ARA*model->VR_prec[k-1];
    #     model->VR_prec[k]=(1-rho_ARA)*model->VR_prec[k-1];
    # } 
    # prov=(1-rho_ARA)*model->A*(model->time[model->k]-model->time[model->k - 1]);
    # //printf("Vright=%f, rho=%f, A=%f, Tk=%f, Tk-1=%f\n",model->Vright,rho,model->A,model->time[model->k],model->time[model->k - 1]);
    # if(nk>0) model->VR_prec[0]=prov;
    # model->Vright+=prov;
    # //printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model->Vright,nk,nk2,model->VR_prec[0],model->VR_prec[1],model->VR_prec[2]); 
    # model->A=pow(rho_QR,f->eval(K)-f->eval(K-1))*model->A;

    # // save old model
    # model->idMod = id;
end