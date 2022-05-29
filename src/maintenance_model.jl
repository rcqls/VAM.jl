abstract type AbstractMaintenanceModel end
init!(mm::AbstractMaintenanceModel) = nothing

mutable struct ARA1 <: AbstractMaintenanceModel
    ρ::Float64
end
params(m::ARA1)::Vector{Float64} = [m.ρ]
params!(m::ARA1, p::Vector{Float64}) = begin;m.ρ = p[1]; nothing; end
nb_params(m::ARA1) = 1
mutable struct ARAInf <: AbstractMaintenanceModel 
	ρ::Float64
end
params(m::ARAInf)::Vector{Float64} = [m.ρ]
params!(m::ARAInf, p::Vector{Float64}) = begin;m.ρ = p[1]; nothing; end
nb_params(m::ARAInf) = 1

mutable struct ARAm <: AbstractMaintenanceModel
    ρ::Float64
    m::Int
end
params(m::ARAm)::Vector{Float64} = [m.ρ]
params!(m::ARAm, p::Vector{Float64}) = begin;m.ρ = p[1]; nothing; end
nb_params(m::ARAm) = 1 # only parameters considered in the optim

struct AGAN <: AbstractMaintenanceModel
end
params(m::AGAN)::Vector{Float64} = []
params!(m::AGAN, p::Vector{Float64}) = nothing
nb_params(m::AGAN) = 0
struct ABAO <: AbstractMaintenanceModel
end
params(m::ABAO)::Vector{Float64} = []
params!(m::ABAO, p::Vector{Float64}) = nothing
nb_params(m::ABAO) = 0

struct AGAP <: AbstractMaintenanceModel
end
params(m::AGAP)::Vector{Float64} = []
params!(m::AGAP, p::Vector{Float64}) = nothing
nb_params(m::AGAP) = 0
mutable struct QAGAN <: AbstractMaintenanceModel
    ρ::Float64
end
params(m::QAGAN)::Vector{Float64} = [m.ρ]
params!(m::QAGAN, p::Vector{Float64}) = begin;m.ρ = p[1]; nothing; end
nb_params(m::QAGAN) = 1

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

abstract type GQRMaintenanceModel <: AbstractMaintenanceModel end
function init!(mm::GQRMaintenanceModel)
    mm.K = 0
end

mutable struct GQR <: GQRMaintenanceModel
    ρ::Float64
    K::Float64
    f::F_GQR
end
params(m::GQR)::Vector{Float64} = [m.ρ]
params!(m::GQR, p::Vector{Float64}) = begin;m.ρ = p[1]; nothing; end
nb_params(m::GQR) = 1

mutable struct GQR_ARA1 <:  GQRMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
end
params(m::GQR_ARA1)::Vector{Float64} = [m.ρQR, m.ρARA]
params!(m::GQR_ARA1, p::Vector{Float64}) = begin; m.ρQR, m.ρARA = p; nothing; end
nb_params(m::GQR_ARA1) = 2
mutable struct GQR_ARAInf <:  GQRMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
end
params(m::GQR_ARAInf)::Vector{Float64} = [m.ρQR, m.ρARA]
params!(m::GQR_ARAInf, p::Vector{Float64}) = begin; m.ρQR, m.ρARA = p; nothing; end
nb_params(m::GQR_ARAInf) = 2
mutable struct GQR_ARAm <: GQRMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
    m::Int
end
params(m::GQR_ARAm)::Vector{Float64} = [m.ρQR, m.ρARA]
params!(m::GQR_ARAm, p::Vector{Float64}) = begin; m.ρQR, m.ρARA = p; nothing; end
nb_params(m::GQR_ARAm) = 2

function update!(m::ARA1, model::AbstractModel; gradient::Bool=false, hessian::Bool=false)
    model.k += 1
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    if hessian
        nb2d = model.nb_params_maintenance * (model.nb_params_maintenance + 1) ÷ 2
        for k in nk:-1:1 
            for i in 0:(model.nb_params_maintenance - 1)
                for j in 0:i
                    model.d2VR_prec[k * nd2d + i * (i + 1) ÷ 2 + j ] = model.d2VR_prec[ (k - 1 ) * nb2d + i * (i + 1) / 2 + j]
                end
            end
        end
        if nk > 0 
            for i in 0:(model.nb_params_maintenance-1)
                for j in 0:i
                    prov = (1-m.ρ)*model.d2A[i * (i + 1) / 2 + j] * (model.time[model.k]-model.time[model.k - 1])
                    model.d2VR_prec[i * (i + 1) / 2 + j] = prov
                    model.d2Vright[i * (i + 1) / 2 + j] += prov
                end
            end
            for j in 0:model.id_params
                prov = model.dA[j] * (model.time[model.k] - model.time[model.k - 1])
                model.d2VR_prec[model.id_params * (model.id_params + 1) / 2 + j] -= prov
                model.d2Vright[model.id_params * (model.id_params + 1) / 2 + j] -= prov
            end
            for i in model.id_params:(model.nb_params_maintenance - 1)
                prov = model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
                model.d2VR_prec[i * (i + 1) / 2 + model.id_params] -= prov
                model.d2Vright[i * (i + 1) / 2 + model.id_params] -= prov
            end
        else
           for i in 0:(model.nb_params_maintenance-1)
                for j in 0:i
                    model.d2Vright[i * (i + 1) / 2 + j] += (1-m.ρ) * model.d2A[i * (i + 1) / 2 + j] * (model.time[model.k]-model.time[model.k - 1])
                end
            end
            for j in 0:model.id_params
                model.d2Vright[model.id_params * (model.id_params + 1) / 2 + j] -= model.dA[j] * (model.time[model.k] - model.time[model.k - 1])
            end
            for i in model.id_params:(model.nb_params_maintenance - 1)
                model.d2Vright[i * (i + 1) / 2 + model.id_params] -= model.dA[i]*(model.time[model.k]-model.time[model.k - 1])
            end 
        end
    end
    if gradient || hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance - 1)
                model.dVR_prec[k * model.nb_params_maintenance + i] = model.dVR_prec[(k-1) * model.nb_params_maintenance + i]
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance - 1)
                prov = (1 - m.ρ) * model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
            prov=model.A * (model.time[model.k] - model.time[model.k - 1])
            model.dVR_prec[model.id_params] -= prov
            model.dVright[model.id_params] -= prov
        else
            for i in 0:(model.nb_params_maintenance - 1)
                model.dVright[i] += (1 - m.ρ) * model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
            end
            model.dVright[model.id_params] -= model.A * (model.time[model.k] - model.time[model.k - 1])
        end
    end
    if nk > 1
        for k in (nk - 1):-1:1
            model.VR_prec[k + 1] = model.VR_prec[k]
        end
    end
    prov = (1 - m.ρ) * model.A * (model.time[model.k] - model.time[model.k - 1])
    if nk > 0
        model.VR_prec[1] = prov
    end
    model.Vright += prov
end

function update!(m::ARAInf, model::AbstractModel; gradient::Bool=false, hessian::Bool=false)
    model.k += 1
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    if hessian
        nb2d = model.nb_params_maintenance * (model.nb_params_maintenance + 1) ÷ 2
        iip = (model.id_params - 1) * model.id_params ÷ 2
        if nk > 1
            for k in nk:-1:1
                for i in 0:(model.nb_params_maintenance - 1)
                    for j in 0:i
                        model.d2VR_prec[k * nb2d + i * (i + 1) ÷ 2+j]=(1 - m.ρ) * model.d2VR_prec[(k - 1) * nb2d + i * (i + 1) ÷ 2 + j]
                    end
                end
                for j in 0:model.id_params
                    model.d2VR_prec[k * nb2d + model.id_params * (model.id_params + 1) ÷ 2 + j] -= model.dVR_prec[(k - 1) * model.nb_params_maintenance + j]
                end
                for i in model.id_params:(model.nb_params_maintenance - 1)
                    model.d2VR_prec[k * nb2d + i * (i + 1) ÷ 2 + model.id_params] -= model.dVR_prec[(k - 1) * model.nb_params_maintenance + i]
                end
            end
        end
        if nk > 1
            for i in 1:model.nb_params_maintenance
                ii = (i - 1) * i ÷ 2
                for j in 1:i
                    model.d2VR_prec[ii + j] = (1-m.ρ) * model.d2A[ii + j] * (model.time[model.k] - model.time[model.k - 1])
                    model.d2Vright[ii + j] = (1 - m.ρ) * model.d2Vleft[ii + j]
                end
            end
            for j in 1:model.id_params
                model.d2VR_prec[iip + j] -= model.dA[j] * (model.time[model.k] - model.time[model.k - 1])
                model.d2Vright[iip + j] -= model.dVleft[j]
            end
            for i in (model.id_params + 1):model.nb_params_maintenance
                ii = (i - 1) * i ÷ 2
                model.d2VR_prec[ii + model.id_params] -= model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
                model.d2Vright[ii + model.id_params] -= model.dVleft[i]
            end
        else
            for i in 1:model.nb_params_maintenance
                ii = (i - 1) * i ÷ 2
                for j in 1:i
                    model.d2Vright[ii + j] = (1 - m.ρ) * model.d2Vleft[ii + j]
                end
            end
            for j in 1:model.id_params
                model.d2Vright[iip + j] -= model.dVleft[j]
            end
            for i in model.id_params:(model.nb_params_maintenance - 1)
                model.d2Vright[(i - 1) * i ÷ 2 + model.id_params] -= model.dVleft[i]
            end
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1 # for(k=nk-1;k>0;k--)
                for i in 1:model.nb_params_maintenance
                    model.dVR_prec[k * model.nb_params_maintenance + i] = (1 - m.ρ) * model.dVR_prec[(k - 1) * model.nb_params_maintenance + i]
                end
                model.dVR_prec[k * model.nb_params_maintenance + model.id_params] -= model.VR_prec[k - 1]
            end
        end
        if nk > 0 
            for i in 1:model.nb_params_maintenance
                model.dVR_prec[i] = (1 - m.ρ) * model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
                model.dVright[i] = (1 - m.ρ) * model.dVleft[i]
            end
            model.dVR_prec[model.id_params] -= model.A * (model.time[model.k] - model.time[model.k - 1])
            model.dVright[model.id_params] -= model.Vleft
        else
            for i in 0:(model.nb_params_maintenance - 1)
                model.dVright[i] = (1-m.ρ) * model.dVleft[i]
            end
            model.dVright[model.id_params] -= model.Vleft
        end
    end
    model.Vright = (1 - m.ρ) * model.Vleft
    if nk > 1
        for k in (nk - 1):-1:1 
            model.VR_prec[k + 1] = (1-m.ρ) * model.VR_prec[k]
        end
    end
    if nk > 0 
        model.VR_prec[1] = (1-m.ρ) * model.A * (model.time[model.k] - model.time[model.k - 1])
    end
end


function update!(m::AGAN, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    model.k += 1
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    model.A = 1
    model.Vright = 0
    for k in 0:(nk - 1)
        model.VR_prec[k + 1] = 0
    end
    if hessian
        for i in 0:(model.nb_params_maintenance - 1)
            model.dA[i] = 0
            model.dVright[i] = 0
            for k in 0:(nk - 1)
                model.dVR_prec[k*model.nb_params_maintenance+i] = 0
            end
            for j in 0:i
                ii = i * (i + 1) ÷ 2
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2A[ii + j] = 0
                model.d2Vright[ii + j] = 0
                for k in 0:(nk - 1)
                    model.d2VR_prec[k * (model.nb_params_maintenance * (model.nb_params_maintenance + 1) ÷ 2) + ii + j] = 0
                end
            end
        end
    end
    if gradient 
        for i in 0:(model.nb_params_maintenance - 1) 
            model.dA[i] = 0
            model.dVright[i] = 0
            for k in 0:(nk - 1)
                model.dVR_prec[k*model.nb_params_maintenance + i] = 0
            end
        end
    end


    # //init QR and GQR type models
    for mm in  model.models
        init!(mm)
    end
end

function update!(m::ABAO, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    if hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance - 1)
                for j in 0:i 
                    jj = (model.nb_params_maintenance * (model.nb_params_maintenance + 1) ÷ 2) + i * (i + 1) ÷ 2 + j
                    model.d2VR_prec[k * jj] = model.d2VR_prec[(k - 1) * jj]
                end
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance - 1)
                ii = i * (i + 1) ÷ 2
                for j in 0:i 
                    prov = model.d2A[ii + j] * (model.time[model.k] - model.time[model.k - 1])
                    model.d2VR_prec[ii + j] = prov
                    model.d2Vright[ii + j] += prov
                end
            end
        else
           for i in 0:(model.nb_params_maintenance - 1) 
                ii = i * (i + 1) ÷ 2
                for j in 0:i 
                    model.d2Vright[ii + j]+=model.d2A[ii + j] * (model.time[model.k] - model.time[model.k - 1])
                end
            end
        end
    end
    if gradient || hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance - 1) 
                model.dVR_prec[k * model.nb_params_maintenance + i] = model.dVR_prec[(k - 1) * model.nb_params_maintenance + i]
            end
        end
        if nk > 0 
            for i in 0:(model.nb_params_maintenance - 1) 
                prov = model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
        else
            for i in 0:(model.nb_params_maintenance - 1) 
                model.dVright[i] += model.dA[i] * (model.time[model.k] - model.time[model.k - 1])
            end
        end
    end
    if nk > 1
        for k in (nk - 1):-1:1
            model.VR_prec[k + 1] = model.VR_prec[k]
        end
    end
    prov = model.A * (model.time[model.k] - model.time[model.k - 1])
    if nk > 0 
        model.VR_prec[1] = prov
    end
    model.Vright += prov
end

function update!(m::AGAP, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
     
    model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    if hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1) 
                 for j in 0:i 
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                end
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1) 
                 for j in 0:i 
                    model.d2VR_prec[i*(i+1)/2+j] = 0
                end
            end
        end
    end
    if gradient || hessian 
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                model.dVR_prec[k*model.nb_params_maintenance+i]=model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1) 
                model.dVR_prec[i]=0;
            end
        end
    end
    if nk > 1
        for k in (n - 1):-1:1
            model.VR_prec[k + 1] = model.VR_prec[k]
        end
    end
    if nk > 0 
        model.VR_prec[0]=0
    end
end


function update!(m::QAGAN, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    model.Vright=0
    for k in 0:(nk - 1)
        model.VR_prec[k + 1]=0
    end
    if hessian
        for i in 0:(model.nb_params_maintenance-1) 
            model.dVright[i] = 0;
            for k in 0:nk
                model.dVR_prec[k*model.nb_params_maintenance+i]=0;
            end
             for j in 0:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2Vright[i*(i+1)/2+j] = 0;
                for k in 0:nk
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=0;
                end
            end
        end
    end
    if gradient 
        for i in 0:(model.nb_params_maintenance-1)
            model.dVright[i] = 0;
            for k in 0:nk
                model.dVR_prec[k*model.nb_params_maintenance+i]=0;
            end
        end
    end
end


function update!(m::QR, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    if  hessian
        for i in 0:(model.nb_params_maintenance-1) 
             for j in 0:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2A[i*(i+1)/2+j] = m.ρ* model.d2A[i*(i+1)/2+j];
                model.d2Vright[i*(i+1)/2+j] = 0;
                for k in 0:nk
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=0;
                end
            end
        end
        for j in 0:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model.d2A[model.id_params*(model.id_params+1)/2+j] = model.d2A[model.id_params*(model.id_params+1)/2+j] + model.dA[j];
        end
        for i in model.id_params:(model.nb_params_maintenance - 1) 
            #//id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model.d2A[i*(i+1)/2+model.id_params] = model.d2A[i*(i+1)/2+model.id_params] + model.dA[i];
        end
    end
    if gradient || hessian 
        for i in 0:(model.nb_params_maintenance - 1)
            model.dA[i] = m.ρ *  model.dA[i]
            model.dVright[i] = 0;
            for k in 0:nk
                model.dVR_prec[k*model.nb_params_maintenance+i]=0
            end
        end
        model.dA[model.id_params] = model.dA[model.id_params] +  model.A;
    end
    model.A=m.ρ * model.A
    model.Vright=0
    for k in 0:(nk - 1)
        model.VR_prec[k + 1]=0
    end
end

function update!(m::GQR, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    model.k += 1;
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    K += 1 
    if hessian
        for i in 0:(model.nb_params_maintenance-1)
             for j in 0:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2A[i*(i+1)/2+j] = pow(m.ρ,f.eval(K)-f.eval(K-1))* model.d2A[i*(i+1)/2+j];
                model.d2Vright[i*(i+1)/2+j] = 0;
                for k in 0:nk
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=0;
                end
            end
        end
        for j in 0:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model.d2A[model.id_params*(model.id_params+1)/2+j] = model.d2A[model.id_params*(model.id_params+1)/2+j] + pow(m.ρ,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ*model.dA[j];
        end
        model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] = model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] + pow(m.ρ,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ*(2*model.dA[model.id_params]+(f.eval(K)-f.eval(K-1)-1)/m.ρ*model.A);
        for i in (model.id_params+1):(model.nb_params_maintenance - 1)
            #//id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model.d2A[i*(i+1)/2+model.id_params] = model.d2A[i*(i+1)/2+model.id_params] + pow(m.ρ,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ*model.dA[i];
        end
    end
    if gradient || hessian
        for i in 0:(model.nb_params_maintenance-1)
            model.dA[i] = pow(m.ρ,f.eval(K)-f.eval(K-1)) *  model.dA[i];
            model.dVright[i] = 0;
            for k in 0:nk
                model.dVR_prec[k*model.nb_params_maintenance+i]=0;
            end
        end
        model.dA[model.id_params] = model.dA[model.id_params] + pow(m.ρ,f.eval(K)-f.eval(K-1)) * (f.eval(K)-f.eval(K-1))/m.ρ * model.A;
    end
    model.A=pow(m.ρ,f.eval(K)-f.eval(K-1))*model.A;
    model.Vright=0;
    for k in 0:nk
        model.VR_prec[k]=0;
    end
end

function update!(m::GQR_ARA1, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    K += 1
    model.k += 1;

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model.k,model.max_mem,nk);
    if hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                end
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    prov= (1-m.ρ_ARA)*model.d2A[i*(i+1)/2+j]*(model.time[model.k]-model.time[model.k - 1]);
                    model.d2VR_prec[i*(i+1)/2+j] = prov;
                    model.d2Vright[i*(i+1)/2+j]+=prov;
                end
            end
            for j in 0:model.id_params
                prov=model.dA[j]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2VR_prec[(model.id_params+1)*(model.id_params+2)/2+j] -= prov;
                model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j] -= prov;
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                prov=model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2VR_prec[i*(i+1)/2+model.id_params+1] -= prov;
                model.d2Vright[i*(i+1)/2+model.id_params+1] -= prov;
            end
        else
           for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2Vright[i*(i+1)/2+j]+=(1-m.ρ_ARA)*model.d2A[i*(i+1)/2+j]*(model.time[model.k]-model.time[model.k - 1]);
                end
            end
            for j in 0:model.id_params
                model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j] -= model.dA[j]*(model.time[model.k]-model.time[model.k - 1]);
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                model.d2Vright[i*(i+1)/2+model.id_params+1] -= model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
            end
        end
        for i in 0:(model.nb_params_maintenance-1)
             for j in 0:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2A[i*(i+1)/2+j] = pow(m.ρ_QR,f.eval(K)-f.eval(K-1))* model.d2A[i*(i+1)/2+j];
            end
        end
        for j in 0:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model.d2A[model.id_params*(model.id_params+1)/2+j] = model.d2A[model.id_params*(model.id_params+1)/2+j] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*model.dA[j];
        end
        model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] = model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*(2*model.dA[model.id_params]+(f.eval(K)-f.eval(K-1)-1)/m.ρ_QR*model.A);
        for i in (model.id_params+1):(model.nb_params_maintenance - 1)
            # //id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model.d2A[i*(i+1)/2+model.id_params] = model.d2A[i*(i+1)/2+model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*model.dA[i];
        end
    end
    if gradient || hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                model.dVR_prec[k*model.nb_params_maintenance+i]=model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                prov=(1-m.ρ_ARA)*model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.dVR_prec[i]=prov;
                model.dVright[i]+=prov;
            end
            prov=model.A*(model.time[model.k]-model.time[model.k - 1]);
            model.dVR_prec[model.id_params+1]-= prov;
            model.dVright[model.id_params+1]-=prov;
        else
            for i in 0:(model.nb_params_maintenance-1)
                model.dVright[i]+=(1-m.ρ_ARA)*model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
            end
            model.dVright[model.id_params+1]-=model.A*(model.time[model.k]-model.time[model.k - 1]);
        end
        for i in 0:(model.nb_params_maintenance-1)
            model.dA[i] = pow(m.ρ_QR,f.eval(K)-f.eval(K-1)) *  model.dA[i];
        end
        model.dA[model.id_params] = model.dA[model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1)) * (f.eval(K)-f.eval(K-1))/m.ρ_QR * model.A;

    end
    for k in nk:-1:1 
        model.VR_prec[k]=model.VR_prec[k-1];
    end
    prov=(1-m.ρ_ARA)*model.A*(model.time[model.k]-model.time[model.k - 1]);
    if nk > 0 
        model.VR_prec[0]=prov
    end
    model.Vright+=prov
    model.A=pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*model.A
end

function update!(m::GQR_ARAInf, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    K += 1
    model.k += 1


    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model.k,model.max_mem,nk);
    if  hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=(1-m.ρ_ARA)*model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                end
            end
            for j in 0:model.id_params
                model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+(model.id_params+1)*(model.id_params+2)/2+j]-= model.dVR_prec[(k-1)*model.nb_params_maintenance+j];
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+model.id_params+1] -= model.dVR_prec[(k-1)*(model.nb_params_maintenance)+i];
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2VR_prec[i*(i+1)/2+j] = (1-m.ρ_ARA)*model.d2A[i*(i+1)/2+j]*(model.time[model.k]-model.time[model.k - 1]);
                    model.d2Vright[i*(i+1)/2+j] = (1-m.ρ_ARA)*model.d2Vleft[i*(i+1)/2+j];
                end
            end
            for j in 0:model.id_params
                model.d2VR_prec[(model.id_params+1)*(model.id_params+2)/2+j] -= model.dA[j]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j] -= model.dVleft[j];
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                model.d2VR_prec[i*(i+1)/2+model.id_params+1] -= model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2Vright[i*(i+1)/2+model.id_params+1] -= model.dVleft[i];
            end
        else
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2Vright[i*(i+1)/2+j] = (1-m.ρ_ARA)*model.d2Vleft[i*(i+1)/2+j];
                end
            end
            for j in 0:model.id_params
                model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j] -= model.dVleft[j];
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                model.d2Vright[i*(i+1)/2+model.id_params+1] -= model.dVleft[i];
            end
        end
        for i in 0:(model.nb_params_maintenance-1)
             for j in 0:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2A[i*(i+1)/2+j] = pow(m.ρ_QR,f.eval(K)-f.eval(K-1))* model.d2A[i*(i+1)/2+j];
            end
        end
        for j in 0:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model.d2A[model.id_params*(model.id_params+1)/2+j] = model.d2A[model.id_params*(model.id_params+1)/2+j] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*model.dA[j];
        end
        model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] = model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*(2*model.dA[model.id_params]+(f.eval(K)-f.eval(K-1)-1)/m.ρ_QR*model.A);
        for i in (model.id_params+1):(model.nb_params_maintenance - 1)
            # //id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model.d2A[i*(i+1)/2+model.id_params] = model.d2A[i*(i+1)/2+model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*model.dA[i];
        end
    end
    if gradient || hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                model.dVR_prec[k*model.nb_params_maintenance+i]=(1-m.ρ_ARA)*model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
            end
            model.dVR_prec[k*model.nb_params_maintenance+model.id_params+1]-=model.VR_prec[k-1];
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                model.dVR_prec[i]=(1-m.ρ_ARA)*model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.dVright[i]=(1-m.ρ_ARA)*model.dVleft[i];
            end
            model.dVR_prec[model.id_params+1]-= model.A*(model.time[model.k]-model.time[model.k - 1]);
            model.dVright[model.id_params+1]-= model.Vleft;
        else
            for i in 0:(model.nb_params_maintenance-1)
                model.dVright[i]=(1-m.ρ_ARA)*model.dVleft[i];
            end
            model.dVright[model.id_params+1]-= model.Vleft;
        end
        for i in 0:(model.nb_params_maintenance-1)
            model.dA[i] = pow(m.ρ_QR,f.eval(K)-f.eval(K-1)) *  model.dA[i];
        end
        model.dA[model.id_params] = model.dA[model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1)) * (f.eval(K)-f.eval(K-1))/m.ρ_QR * model.A;
    end
    model.Vright=(1-m.ρ_ARA)*model.Vleft;
    for k in nk:-1:1 
        model.VR_prec[k]=(1-m.ρ_ARA)*model.VR_prec[k-1]
    end
    if nk > 0 
        model.VR_prec[0]=(1-m.ρ_ARA)*model.A*(model.time[model.k]-model.time[model.k - 1])
    end
    model.A=pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*model.A
end

function update!(m::ARAm, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    model.k += 1;

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    nk2 = nk
    if nk > m-1 
        nk2 = m-1
    end

    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model.k,model.max_mem,nk);
    if hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                end
            end
        end

        if  (model.k >= m) && (nk2 > 0)
            if nk > nk2
                for i in 0:(model.nb_params_maintenance-1)
                     for j in 0:i
                        model.d2Vright[i*(i+1)/2+j]-=m.ρ*model.d2VR_prec[(nk2-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                        model.d2VR_prec[nk2*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=(1-m.ρ)*model.d2VR_prec[(nk2-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                    end
                end
                for j in 0:model.id_params
                    model.d2Vright[model.id_params*(model.id_params+1)/2+j]-= model.dVR_prec[(nk2-1)*model.nb_params_maintenance+j];
                    model.d2VR_prec[nk2*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+model.id_params*(model.id_params+1)/2+j]-= model.dVR_prec[(nk2-1)*model.nb_params_maintenance+j];
                end
                for i in model.id_params:(model.nb_params_maintenance - 1) 
                    model.d2Vright[i*(i+1)/2+model.id_params] -= model.dVR_prec[(nk2-1)*(model.nb_params_maintenance)+i];
                    model.d2VR_prec[nk2*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+model.id_params] -= model.dVR_prec[(nk2-1)*(model.nb_params_maintenance)+i];
                end
            else
                for i in 0:(model.nb_params_maintenance-1)
                     for j in 0:i 
                        model.d2Vright[i*(i+1)/2+j]-=m.ρ*model.d2VR_prec[(nk2-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]
                    end
                end
                for j in 0:model.id_params 
                    model.d2Vright[model.id_params*(model.id_params+1)/2+j]-= model.dVR_prec[(nk2-1)*model.nb_params_maintenance+j];
                end
                for i in model.id_params:(model.nb_params_maintenance - 1)  
                    model.d2Vright[i*(i+1)/2+model.id_params] -= model.dVR_prec[(nk2-1)*(model.nb_params_maintenance)+i];
                end
            end
        end

        for k in nk2:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2Vright[i*(i+1)/2+j]-=m.ρ*model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=(1-m.ρ)*model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                end
            end
            for j in 0:model.id_params
                model.d2Vright[model.id_params*(model.id_params+1)/2+j]-= model.dVR_prec[(k-1)*model.nb_params_maintenance+j];
                model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+model.id_params*(model.id_params+1)/2+j]-= model.dVR_prec[(k-1)*model.nb_params_maintenance+j];
            end
            for i in model.id_params:(model.nb_params_maintenance - 1) 
                model.d2Vright[i*(i+1)/2+model.id_params] -= model.dVR_prec[(k-1)*(model.nb_params_maintenance)+i];
                model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+model.id_params] -= model.dVR_prec[(k-1)*(model.nb_params_maintenance)+i];
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    prov= (1-m.ρ)*model.d2A[i*(i+1)/2+j]*(model.time[model.k]-model.time[model.k - 1]);
                    model.d2VR_prec[i*(i+1)/2+j] = prov;
                    model.d2Vright[i*(i+1)/2+j]+=prov;
                end
            end
            for j in 0:model.id_params
                prov=model.dA[j]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2VR_prec[model.id_params*(model.id_params+1)/2+j] -= prov;
                model.d2Vright[model.id_params*(model.id_params+1)/2+j] -= prov;
            end
            for i in model.id_params:(model.nb_params_maintenance - 1) 
                prov=model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2VR_prec[i*(i+1)/2+model.id_params] -= prov;
                model.d2Vright[i*(i+1)/2+model.id_params] -= prov;
            end
        else
           for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2Vright[i*(i+1)/2+j]+=(1-m.ρ)*model.d2A[i*(i+1)/2+j]*(model.time[model.k]-model.time[model.k - 1]);
                end
            end
            for j in 0:model.id_params
                model.d2Vright[model.id_params*(model.id_params+1)/2+j] -= model.dA[j]*(model.time[model.k]-model.time[model.k - 1]);
            end
            for i in model.id_params:(model.nb_params_maintenance - 1) 
                model.d2Vright[i*(i+1)/2+model.id_params] -= model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
            end
        end
    end
    if gradient || hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                model.dVR_prec[k*model.nb_params_maintenance+i]=model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
            end
        end
        if (model.k >= m) && (nk2 > 0) 
            if nk > nk2
                for i in 0:(model.nb_params_maintenance-1)
                    model.dVright[i]-=m.ρ*model.dVR_prec[(nk2-1)*model.nb_params_maintenance+i];
                    model.dVR_prec[nk2*model.nb_params_maintenance+i]=(1-m.ρ)*model.dVR_prec[(nk2-1)*model.nb_params_maintenance+i];
                end
                model.dVR_prec[nk2 * (model.nb_params_maintenance)+model.id_params] -= model.VR_prec[nk2-1];
                model.dVright[model.id_params]-=model.VR_prec[nk2-1];
            else
                for i in 0:(model.nb_params_maintenance-1) 
                    model.dVright[i]-=m.ρ*model.dVR_prec[(nk2-1)*model.nb_params_maintenance+i];
                end
                model.dVright[model.id_params]-=model.VR_prec[nk2-1];
            end
        end
        for k in nk2:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                model.dVright[i]-=m.ρ*model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
                model.dVR_prec[k*model.nb_params_maintenance+i]=(1-m.ρ)*model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
            end
            model.dVR_prec[k * (model.nb_params_maintenance)+model.id_params] -= model.VR_prec[k-1];
            model.dVright[model.id_params]-=model.VR_prec[k-1];
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                prov=(1-m.ρ)*model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.dVR_prec[i]=prov;
                model.dVright[i]+=prov;
            end
            prov=model.A*(model.time[model.k]-model.time[model.k - 1]);
            model.dVR_prec[model.id_params]-= prov;
            model.dVright[model.id_params]-=prov;
        else
            for i in 0:(model.nb_params_maintenance-1)
                model.dVright[i]+=(1-m.ρ)*model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
            end
            model.dVright[model.id_params]-=model.A*(model.time[model.k]-model.time[model.k - 1]);
        end
    end
    # //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
    for k=nk:-1:nk2
        model.VR_prec[k]=model.VR_prec[k-1];
    end
    if (model.k >= m) && (nk2 > 0)
        if nk > nk2
            model.Vright-=m.ρ*model.VR_prec[nk2-1];
            model.VR_prec[nk2]=(1-m.ρ)*model.VR_prec[nk2-1];
        else 
            model.Vright-=m.ρ*model.VR_prec[nk2-1];
        end
    end
    for k=nk2:-1:1
        model.Vright-=m.ρ*model.VR_prec[k-1];
        model.VR_prec[k]=(1-m.ρ)*model.VR_prec[k-1];
    end
    prov=(1-m.ρ)*model.A*(model.time[model.k]-model.time[model.k - 1]);
    #//printf("Vright=%f, m.ρ=%f, A=%f, Tk=%f, Tk-1=%f\n",model.Vright,m.ρ,model.A,model.time[model.k],model.time[model.k - 1]);
    if nk > 0 
        model.VR_prec[0]=prov
    end
    model.Vright+=prov;
    #//printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model.Vright,nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
end

function update!(m::GQR_ARAm, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    model.k += 1;
    K += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    nk2 = nk
    if nk > m-1
        nk2 = m-1
    end
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model.k,model.max_mem,nk);
    if hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                end
            end
        end

        if  (model.k >= m) && (nk2 > 0)
            if nk > nk2
                for i in 0:(model.nb_params_maintenance-1)
                     for j in 0:i
                        model.d2Vright[i*(i+1)/2+j]-=m.ρ_ARA*model.d2VR_prec[(nk2-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                        model.d2VR_prec[nk2*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=(1-m.ρ_ARA)*model.d2VR_prec[(nk2-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                    end
                end
                for j in 0:model.id_params
                    model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j]-= model.dVR_prec[(nk2-1)*model.nb_params_maintenance+j];
                    model.d2VR_prec[nk2*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+(model.id_params+1)*(model.id_params+2)/2+j]-= model.dVR_prec[(nk2-1)*model.nb_params_maintenance+j];
                end
                for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                    model.d2Vright[i*(i+1)/2+model.id_params+1] -= model.dVR_prec[(nk2-1)*(model.nb_params_maintenance)+i];
                    model.d2VR_prec[nk2*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+model.id_params+1] -= model.dVR_prec[(nk2-1)*(model.nb_params_maintenance)+i];
                end
            else
                for i in 0:(model.nb_params_maintenance-1)
                     for j in 0:i 
                        model.d2Vright[i*(i+1)/2+j]-=m.ρ_ARA*model.d2VR_prec[(nk2-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                    end
                end
                for j in 0:(model.id_params+1) 
                    model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j]-= model.dVR_prec[(nk2-1)*model.nb_params_maintenance+j]
                end
                for i in (model.id_params+1):(model.nb_params_maintenance - 1) 
                    model.d2Vright[i*(i+1)/2+model.id_params+1] -= model.dVR_prec[(nk2-1)*(model.nb_params_maintenance)+i];
                end
            end
        end

        for k in nk2:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2Vright[i*(i+1)/2+j]-=m.ρ_ARA*model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                    model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j]=(1-m.ρ_ARA)*model.d2VR_prec[(k-1)*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+j];
                end
            end
            for j in 0:model.id_params
                model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j]-= model.dVR_prec[(k-1)*model.nb_params_maintenance+j];
                model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+(model.id_params+1)*(model.id_params+2)/2+j]-= model.dVR_prec[(k-1)*model.nb_params_maintenance+j];
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                model.d2Vright[i*(i+1)/2+model.id_params+1] -= model.dVR_prec[(k-1)*(model.nb_params_maintenance)+i];
                model.d2VR_prec[k*(model.nb_params_maintenance*(model.nb_params_maintenance+1)/2)+i*(i+1)/2+model.id_params+1] -= model.dVR_prec[(k-1)*(model.nb_params_maintenance)+i];
            end
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    prov= (1-m.ρ_ARA)*model.d2A[i*(i+1)/2+j]*(model.time[model.k]-model.time[model.k - 1]);
                    model.d2VR_prec[i*(i+1)/2+j] = prov;
                    model.d2Vright[i*(i+1)/2+j]+=prov;
                end
            end
            for j in 0:model.id_params
                prov=model.dA[j]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2VR_prec[(model.id_params+1)*(model.id_params+2)/2+j] -= prov;
                model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j] -= prov;
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                prov=model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.d2VR_prec[i*(i+1)/2+model.id_params+1] -= prov;
                model.d2Vright[i*(i+1)/2+model.id_params+1] -= prov;
            end
        else
           for i in 0:(model.nb_params_maintenance-1)
                 for j in 0:i
                    model.d2Vright[i*(i+1)/2+j]+=(1-m.ρ_ARA)*model.d2A[i*(i+1)/2+j]*(model.time[model.k]-model.time[model.k - 1]);
                end
            end
            for j in 0:model.id_params
                model.d2Vright[(model.id_params+1)*(model.id_params+2)/2+j] -= model.dA[j]*(model.time[model.k]-model.time[model.k - 1]);
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                model.d2Vright[i*(i+1)/2+model.id_params+1] -= model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
            end
        end
        for i in 0:(model.nb_params_maintenance-1)
             for j in 0:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2A[i*(i+1)/2+j] = pow(m.ρ_QR,f.eval(K)-f.eval(K-1))* model.d2A[i*(i+1)/2+j];
            end
        end
        for j in 0:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model.d2A[model.id_params*(model.id_params+1)/2+j] = model.d2A[model.id_params*(model.id_params+1)/2+j] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*model.dA[j];
        end
        model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] = model.d2A[model.id_params*(model.id_params+1)/2+model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*(2*model.dA[model.id_params]+(f.eval(K)-f.eval(K-1)-1)/m.ρ_QR*model.A);
        for i in (model.id_params+1):(model.nb_params_maintenance - 1)
            # //id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model.d2A[i*(i+1)/2+model.id_params] = model.d2A[i*(i+1)/2+model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*(f.eval(K)-f.eval(K-1))/m.ρ_QR*model.dA[i];
        end
    end
    if gradient || hessian
        for k in nk:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                model.dVR_prec[k*model.nb_params_maintenance+i]=model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
            end
        end
        if  (model.k >= m) && (nk2 > 0)
            if nk > nk2
                for i in 0:(model.nb_params_maintenance-1)
                    model.dVright[i]-=m.ρ_ARA*model.dVR_prec[(nk2-1)*model.nb_params_maintenance+i];
                    model.dVR_prec[nk2*model.nb_params_maintenance+i]=(1-m.ρ_ARA)*model.dVR_prec[(nk2-1)*model.nb_params_maintenance+i];
                end
                model.dVR_prec[nk2 * (model.nb_params_maintenance)+model.id_params+1] -= model.VR_prec[nk2-1];
                model.dVright[model.id_params+1]-=model.VR_prec[nk2-1];
            else
                for i in 0:(model.nb_params_maintenance-1) 
                    model.dVright[i]-=m.ρ_ARA*model.dVR_prec[(nk2-1)*model.nb_params_maintenance+i]
                end
                model.dVright[model.id_params+1]-=model.VR_prec[nk2-1];
            end
        end
        for k in nk2:-1:1
            for i in 0:(model.nb_params_maintenance-1)
                model.dVright[i]-=m.ρ_ARA*model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
                model.dVR_prec[k*model.nb_params_maintenance+i]=(1-m.ρ_ARA)*model.dVR_prec[(k-1)*model.nb_params_maintenance+i];
            end
            model.dVR_prec[k * (model.nb_params_maintenance)+model.id_params+1] -= model.VR_prec[k-1];
            model.dVright[model.id_params+1]-=model.VR_prec[k-1];
        end
        if nk > 0
            for i in 0:(model.nb_params_maintenance-1)
                prov=(1-m.ρ_ARA)*model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
                model.dVR_prec[i]=prov;
                model.dVright[i]+=prov;
            end
            prov=model.A*(model.time[model.k]-model.time[model.k - 1]);
            model.dVR_prec[model.id_params+1]-= prov;
            model.dVright[model.id_params+1]-=prov;
        else
            for i in 0:(model.nb_params_maintenance-1)
                model.dVright[i]+=(1-m.ρ_ARA)*model.dA[i]*(model.time[model.k]-model.time[model.k - 1]);
            end
            model.dVright[model.id_params+1]-=model.A*(model.time[model.k]-model.time[model.k - 1]);
        end
        for i in 0:(model.nb_params_maintenance-1)
            model.dA[i] = pow(m.ρ_QR,f.eval(K)-f.eval(K-1)) *  model.dA[i];
        end
        model.dA[model.id_params] = model.dA[model.id_params] + pow(m.ρ_QR,f.eval(K)-f.eval(K-1)) * (f.eval(K)-f.eval(K-1))/m.ρ_QR * model.A;
    end
    # //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
    for k=nk:-1:nk2
        model.VR_prec[k]=model.VR_prec[k-1];
    end
    if (model.k >= m) && (nk2 > 0)
        if nk > nk2
            model.Vright-=m.ρ_ARA*model.VR_prec[nk2-1];
            model.VR_prec[nk2]=(1-m.ρ_ARA)*model.VR_prec[nk2-1];
        else 
            model.Vright-=m.ρ_ARA*model.VR_prec[nk2-1]
        end
    end
    for k=nk2:-1:1
        model.Vright-=m.ρ_ARA*model.VR_prec[k-1];
        model.VR_prec[k]=(1-m.ρ_ARA)*model.VR_prec[k-1];
    end
    prov=(1-m.ρ_ARA)*model.A*(model.time[model.k]-model.time[model.k - 1]);
    #//printf("Vright=%f, m.ρ=%f, A=%f, Tk=%f, Tk-1=%f\n",model.Vright,m.ρ,model.A,model.time[model.k],model.time[model.k - 1]);
    if nk > 0 
        model.VR_prec[0]=prov
    end
    model.Vright+=prov;
    #//printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model.Vright,nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
    model.A=pow(m.ρ_QR,f.eval(K)-f.eval(K-1))*model.A;
end