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

abstract type GQRMaintenanceModel <: AbstractMaintenanceModel end
function init!(mm::GQRMaintenanceModel)
    mm.K = 0
end

mutable struct GQR <: GQRMaintenanceModel
    ρ::Float64
    f::Function
    K::Float64
end
GQR(ρ::Float64, f::Function) = GQR(ρ, f, 0)
params(m::GQR)::Vector{Float64} = [m.ρ]
params!(m::GQR, p::Vector{Float64}) = begin;m.ρ = p[1]; nothing; end
nb_params(m::GQR) = 1

mutable struct GQR_ARA1 <:  GQRMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    f::Function
    K::Float64
end
GQR_ARA1(ρQR::Float64, ρARA::Float64, f::Function) = GQR_ARA1(ρQR, ρARA, f, 0)
params(m::GQR_ARA1)::Vector{Float64} = [m.ρQR, m.ρARA]
params!(m::GQR_ARA1, p::Vector{Float64}) = begin; m.ρQR, m.ρARA = p; nothing; end
nb_params(m::GQR_ARA1) = 2
mutable struct GQR_ARAInf <:  GQRMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    f::Function
    K::Float64
end
GQR_ARAInf(ρQR::Float64, ρARA::Float64, f::Function) = GQR_ARAInf(ρQR, ρARA, f, 0)
params(m::GQR_ARAInf)::Vector{Float64} = [m.ρQR, m.ρARA]
params!(m::GQR_ARAInf, p::Vector{Float64}) = begin; m.ρQR, m.ρARA = p; nothing; end
nb_params(m::GQR_ARAInf) = 2
mutable struct GQR_ARAm <: GQRMaintenanceModel
    ρQR::Float64
    ρARA::Float64
    m::Int
    f::Function
    K::Float64
end
GQR_ARAm(ρQR::Float64, ρARA::Float64, m::Int, f::Function) = GQR_ARAm(ρQR, ρARA, m, f, 0)
params(m::GQR_ARAm)::Vector{Float64} = [m.ρQR, m.ρARA]
params!(m::GQR_ARAm, p::Vector{Float64}) = begin; m.ρQR, m.ρARA = p; nothing; end
nb_params(m::GQR_ARAm) = 2

function update!(m::ARA1, model::AbstractModel; gradient::Bool=false, hessian::Bool=false)
    inc!(model) #model.k += 1
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    npm = model.nb_params_maintenance
    if hessian
        npm2 = ind_nb(npm)
        if nk > 1
            for k in (nk - 1):-1:1 
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2VR_prec[k * nd2d + ij ] = model.d2VR_prec[ (k - 1 ) * npm2 + ij]
                    end
                end
            end
        end
        if nk > 0 
            for i in 1:npm
                for j in 1:i
                    ij = ind_ij(i, j)
                    prov = (1-m.ρ) * model.d2A[ij] * model.Δt
                    model.d2VR_prec[ij] = prov
                    model.d2Vright[ij] += prov
                end
            end
            for j in 1:model.id_params
                jj = ind_ij(model.id_params, j)
                prov = model.dA[j] * model.Δt
                model.d2VR_prec[jj] -= prov
                model.d2Vright[jj] -= prov
            end
            for i in (model.id_params + 1):npm
                ii = ind_ij(i, model.id_params)
                prov = model.dA[i] * model.Δt
                model.d2VR_prec[ii] -= prov
                model.d2Vright[ii] -= prov
            end
        else
           for i in 1:npm
                for j in 1:i
                    model.d2Vright[ind_ij(i, j)] += (1-m.ρ) * model.d2A[ind_ij(i, j)] * model.Δt
                end
            end
            for j in 1:model.id_params
                model.d2Vright[ind_ij(model.id_params, j)] -= model.dA[j] * model.Δt
            end
            for i in (model.id_params + 1):model.nb_params_maintenance
                model.d2Vright[ind_ij(i, model.id_params)] -= model.dA[i]*model.Δt
            end 
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    model.dVR_prec[k * npm + i] = model.dVR_prec[(k-1) * npm + i]
                end
            end
        end
        if nk > 0
            for i in 1:npm
                prov = (1 - m.ρ) * model.dA[i] * model.Δt
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
            prov=model.A * model.Δt
            model.dVR_prec[model.id_params] -= prov
            model.dVright[model.id_params] -= prov
        else
            for i in 1:npm
                model.dVright[i] += (1 - m.ρ) * model.dA[i] * model.Δt
            end
            model.dVright[model.id_params] -= model.A * model.Δt
        end
    end
    if nk > 1
        for k in (nk - 1):-1:1
            model.VR_prec[k + 1] = model.VR_prec[k]
        end
    end
    prov = (1 - m.ρ) * model.A * model.Δt
    if nk > 0
        model.VR_prec[1] = prov
    end
    model.Vright += prov
end

function update!(m::ARAInf, model::AbstractModel; gradient::Bool=false, hessian::Bool=false)
    inc!(model) #model.k += 1
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    npm = model.nb_params_maintenance
    if hessian
        npm2 = ind_nb(npm)
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2VR_prec[k * npm2 + ij] = (1 - m.ρ) * model.d2VR_prec[(k - 1) * npm2 + ij]
                    end
                end
                for j in 1:model.id_params
                    model.d2VR_prec[k * npm2 + ind_ij(model.id_params, j)] -= model.dVR_prec[(k - 1) * npm + j]
                end
                for i in (model.id_params + 1):npm
                    model.d2VR_prec[k * npm2 + ind_ij(i, model.id_params)] -= model.dVR_prec[(k - 1) * npm + i]
                end
            end
        end
        if nk > 0
            for i in 1:npm
                for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2VR_prec[ij] = (1 - m.ρ) * model.d2A[ij] * model.Δt
                    model.d2Vright[ij] = (1 - m.ρ) * model.d2Vleft[ij]
                end
            end
            for j in 1:model.id_params
                jj = ind_ij(model.id_params, j)
                model.d2VR_prec[jj] -= model.dA[j] * model.Δt
                model.d2Vright[jj] -= model.dVleft[j]
            end
            for i in (model.id_params + 1):npm
                ii = ind_ij(i, model.id_params)
                model.d2VR_prec[ii] -= model.dA[i] * model.Δt
                model.d2Vright[ii] -= model.dVleft[i]
            end
        else
            for i in 1:npm
                for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2Vright[ij] = (1 - m.ρ) * model.d2Vleft[ij]
                end
            end
            for j in 1:model.id_params
                jj = ind_ij(model.id_params, j)
                model.d2Vright[jj] -= model.dVleft[j]
            end
            for i in (model.id_params + 1):npm
                model.d2Vright[ind_ij(i, model.id_params)] -= model.dVleft[i]
            end
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    model.dVR_prec[k * npm + i] = (1 - m.ρ) * model.dVR_prec[(k - 1) * npm + i]
                end
                model.dVR_prec[k * npm + model.id_params] -= model.VR_prec[k - 1]
            end
        end
        if nk > 0 
            for i in 1:npm
                model.dVR_prec[i] = (1 - m.ρ) * model.dA[i] * model.Δt
                model.dVright[i] = (1 - m.ρ) * model.dVleft[i]
            end
            model.dVR_prec[model.id_params] -= model.A * model.Δt
            model.dVright[model.id_params] -= model.Vleft
        else
            for i in 1:npm
                model.dVright[i] = (1 - m.ρ) * model.dVleft[i]
            end
            model.dVright[model.id_params] -= model.Vleft
        end
    end
    model.Vright = (1 - m.ρ) * model.Vleft
    if nk > 1
        for k in (nk - 1):-1:1 
            model.VR_prec[k + 1] = (1 - m.ρ) * model.VR_prec[k]
        end
    end
    if nk > 0 
        model.VR_prec[1] = (1 - m.ρ) * model.A * model.Δt
    end
end


function update!(m::AGAN, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    inc!(model) #model.k += 1
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    model.A = 1
    model.Vright = 0
    fill!(model.VR_prec, 0.0)
    if hessian || gradient
        ## short replacement
        fill!(model.dA, 0.0)
        fill!(model.dVright, 0.0)
        fill!(model.dVR_prec, 0.0)
    end
    if hessian
        fill!(model.d2A, 0.0)
        fill!(model.d2Vright, 0.0)
        fill!(model.d2VR_prec, 0.0)
    end
    # //init QR and GQR type models
    for mm in  model.models
        init!(mm)
    end
end

function update!(m::ABAO, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    inc!(model) #model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    npm = model.nb_params_maintenance
    if hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    for j in 1:i 
                        ij = ind_ij(npm, 0) + ind_ij(i, j)
                        model.d2VR_prec[k * ij] = model.d2VR_prec[(k - 1) * ij]
                    end
                end
            end
        end
        if nk > 0
            for i in 1:npm
                for j in 1:i
                    ij = ind_ij(i, j)
                    prov = model.d2A[ij] * model.Δt
                    model.d2VR_prec[ij] = prov
                    model.d2Vright[ij] += prov
                end
            end
        else
           for i in 1:npm
                for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2Vright[ij] += model.d2A[ij] * model.Δt
                end
            end
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm 
                    model.dVR_prec[k * npm + i] = model.dVR_prec[(k - 1) * npm + i]
                end
            end
        end
        if nk > 0 
            for i in 1:npm 
                prov = model.dA[i] * model.Δt
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
        else
            for i in 1:npm
                model.dVright[i] += model.dA[i] * model.Δt
            end
        end
    end
    if nk > 1
        for k in (nk - 1):-1:1
            model.VR_prec[k + 1] = model.VR_prec[k]
        end
    end
    prov = model.A * model.Δt
    if nk > 0 
        model.VR_prec[1] = prov
    end
    model.Vright += prov
end

function update!(m::AGAP, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
     
    inc!(model) #model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    npm = model.nb_params_maintenance
    if hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    for j in 1:i
                        model.d2VR_prec[k * ind_nb(npm) + ind_ij(i, j)] = model.d2VR_prec[(k-1) * ind_nb(npm) + ind_ij(i, j)]
                    end
                end
            end
        end
        if nk > 0
            for i in 1:npm
                for j in 1:i
                    model.d2VR_prec[ind_ij(i, j)] = 0
                end
            end
        end
    end
    if gradient || hessian 
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    model.dVR_prec[k * npm + i] = model.dVR_prec[(k - 1) * npm + i]
                end
            end
        end
        if nk > 0
            for i in 1:npm 
                model.dVR_prec[i] = 0
            end
        end
    end
    if nk > 1
        for k in (n - 1):-1:1
            model.VR_prec[k + 1] = model.VR_prec[k]
        end
    end
    if nk > 0 
        model.VR_prec[0] = 0
    end
end


function update!(m::QAGAN, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    inc!(model) #model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    model.Vright=0
    for k in 0:(nk - 1)
        model.VR_prec[k + 1]=0
    end
    npm = model.nb_params_maintenance
    if hessian
        for i in 1:npm 
            model.dVright[i] = 0
            for k in 1:nk
                model.dVR_prec[k * npm + i] = 0
            end
             for j in 1:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                ij = ind_ij(i, j)
                model.d2Vright[ij] = 0
                for k in 1:nk
                    model.d2VR_prec[k * ind_nb(npm) + ij] = 0
                end
            end
        end
    end
    if gradient 
        for i in 1:npm
            model.dVright[i] = 0
            for k in 1:nk
                model.dVR_prec[k * npm + i] = 0
            end
        end
    end
end


function update!(m::QR, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    inc!(model) #model.k += 1

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    npm = model.nb_params_maintenance
    if hessian
        for i in 1:npm 
             for j in 1:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                ij = ind_ij(i, j)
                model.d2A[ij] = m.ρ * model.d2A[ij]
                model.d2Vright[ij] = 0
                for k in 1:nk
                    model.d2VR_prec[k * ind_nb(npm) + ij] = 0
                end
            end
        end
        for j in 1:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            jj = ind_ij(model.id_params, j)
            model.d2A[jj] = model.d2A[jj] + model.dA[j]
        end
        for i in (model.id_params + 1):npm
            #//id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ii = ind_ij(i, model.id_params)
            model.d2A[ii] = model.d2A[ii] + model.dA[i]
        end
    end
    if gradient || hessian 
        for i in 1:npm
            model.dA[i] = m.ρ *  model.dA[i]
            model.dVright[i] = 0
            for k in 1:nk
                model.dVR_prec[k * npm + i] = 0
            end
        end
        model.dA[model.id_params] = model.dA[model.id_params] +  model.A
    end
    model.A = m.ρ * model.A
    model.Vright=0
    for k in 1:nk
        model.VR_prec[k] = 0
    end
end

function update!(m::GQR, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    inc!(model) #model.k += 1
    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    K += 1
    npm = model.nb_params_maintenance
    δ = (m.f(K) - m.f(K-1))
    if hessian
        for i in 1:npm
             for j in 1:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                ij = ind_ij(i, j)
                model.d2A[ij] = m.ρ^δ  * model.d2A[ij]
                model.d2Vright[ij] = 0
                for k in 1:nk
                    model.d2VR_prec[k * ind_nb(npm) + ij] = 0
                end
            end
        end
        for j in 1:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            jj = ind_ij(model.id_params, j)
            model.d2A[jj] = model.d2A[jj] + m.ρ^δ * δ/m.ρ * model.dA[j]
        end
        model.d2A[ind_ij(model.id_params, model.id_params)] += m.ρ^δ * δ/m.ρ * (2 * model.dA[model.id_params] + (δ - 1)/m.ρ * model.A)
        for i in (model.id_params+1):npm
            #//id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ii = ind_ij(i, model.id_params)
            model.d2A[ii] += m.ρ^δ * δ / m.ρ * model.dA[i]
        end
    end
    if gradient || hessian
        for i in 1:npm
            model.dA[i] = m.ρ^δ * model.dA[i]
            model.dVright[i] = 0
            fill!(model.dVR_prec, 0.0)
        end
        model.dA[model.id_params] = model.dA[model.id_params] + m.ρ^δ * δ / m.ρ * model.A
    end
    model.A = m.ρ^δ * model.A
    model.Vright=0
    for k in 0:nk
        model.VR_prec[k]=0
    end
end

function update!(m::GQR_ARA1, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    K += 1
    inc!(model) #model.k += 1;

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model.k,model.max_mem,nk);
    npm = model.nb_params_maintenance
    δ = (m.f(K) - m.f(K-1))
    if hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2VR_prec[k * ind_nb(npm) + ij] = model.d2VR_prec[(k - 1) * ind_nb(npm) + ij]
                    end
                end
            end
        end
        if nk > 0
            for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    prov= (1 - m.ρ_ARA) * model.d2A[ij] * model.Δt
                    model.d2VR_prec[ij] = prov
                    model.d2Vright[ij] += prov
                end
            end
            for j in 1:model.id_params
                prov = model.dA[j] * model.Δt
                jj = ind_ij(model.id_params + 1, j)
                model.d2VR_prec[jj] -= prov
                model.d2Vright[jj] -= prov
            end
            for i in (model.id_params + 2):npm
                ii = ind_ij(i, model.id_params + 1)
                prov = model.dA[i] * model.Δt
                model.d2VR_prec[ii] -= prov
                model.d2Vright[ii] -= prov
            end
        else
           for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2Vright[ij] += (1 - m.ρ_ARA) * model.d2A[ij] * model.Δt
                end
            end
            for j in 1:model.id_params
                jj = ind_ij(model.id_params + 1, j)
                model.d2Vright[jj] -= model.dA[j] * model.Δt
            end
            for i in (model.id_params+1):npm
                ii = ind_ij(i, model.id_params + 1)
                model.d2Vright[ii] -= model.dA[i] * model.Δt
            end
        end
        for i in 1:npm
             for j in 1:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                ij = ind_ij(i, j)
                model.d2A[ij] = m.ρ_QR^δ * model.d2A[ij]
            end
        end
        for j in 1:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            jj = ind_ij(model.id_params, j)
            model.d2A[jj] += m.ρ_QR^δ * δ / m.ρ_QR * model.dA[j]
        end
        model.d2A[ind_ij(model.id_params, model.id_params)] += m.ρ_QR^δ * δ / m.ρ_QR * (2 * model.dA[model.id_params] + (δ - 1) / m.ρ_QR * model.A)
        for i in (model.id_params+1):npm
            # //id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ii = ind_ij(i, model.id_params)
            model.d2A[ii] += m.ρ_QR^δ * δ / m.ρ_QR * model.dA[i]
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    model.dVR_prec[k * npm + i] = model.dVR_prec[(k-1) * npm + i]
                end
            end
        end
        if nk > 0
            for i in 1:npm
                prov = (1-m.ρ_ARA) * model.dA[i] * model.Δt
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
            prov = model.A * model.Δt
            model.dVR_prec[model.id_params + 1] -= prov
            model.dVright[model.id_params + 1] -= prov
        else
            for i in 1:npm
                model.dVright[i] += (1 - m.ρ_ARA) * model.dA[i] * model.Δt
            end
            model.dVright[model.id_params + 1] -= model.A * model.Δt
        end
        for i in 1:npm
            model.dA[i] = m.ρ_QR^δ * model.dA[i]
        end
        model.dA[model.id_params] += m.ρ_QR^δ * δ / m.ρ_QR * model.A

    end
    if nk > 1
        for k in (nk - 1):-1:1 
            model.VR_prec[k + 1] = model.VR_prec[k]
        end
    end
    prov = (1 - m.ρ_ARA) * model.A * model.Δt
    if nk > 0 
        model.VR_prec[1] = prov
    end
    model.Vright += prov
    model.A *= m.ρ_QR^δ
end

function update!(m::GQR_ARAInf, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    K += 1
    inc!(model) #model.k += 1


    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model.k,model.max_mem,nk);
    npm = model.nb_params_maintenance
    δ = (m.f(K) - m.f(K-1))
    if hessian
        npm2 = ind_nb(npm)
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2VR_prec[k * npm2 + ij] = (1 - m.ρ_ARA) * model.d2VR_prec[(k - 1) * npm2 + ij]
                    end
                end
            end
            for j in 1:model.id_params
                model.d2VR_prec[k * npm2 + ind_ij(model.id_params+1, j)] -= model.dVR_prec[(k-1) * npm + j]
            end
            for i in (model.id_params + 2):npm
                model.d2VR_prec[k * npm2 + ind_ij(i, model.id_params+1)] -= model.dVR_prec[(k-1) * npm + i]
            end
        end
        if nk > 0
            for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2VR_prec[ij] = (1 - m.ρ_ARA) * model.d2A[ij] * model.Δt
                    model.d2Vright[ij] = (1 - m.ρ_ARA) * model.d2Vleft[ij]
                end
            end
            for j in 1:model.id_params
                jj = ind_ij(model.id_params + 1, j)
                model.d2VR_prec[jj] -= model.dA[j] * model.Δt
                model.d2Vright[jj] -= model.dVleft[j]
            end
            for i in (model.id_params+1):npm
                ii = ind_ij(i, model.id_params + 1)
                model.d2VR_prec[ii] -= model.dA[i] * model.Δt
                model.d2Vright[ii] -= model.dVleft[i]
            end
        else
            for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2Vright[ij] = (1-m.ρ_ARA) * model.d2Vleft[ij]
                end
            end
            for j in 1:model.id_params
                jj = ind_ij(model.id_params + 1, j)
                model.d2Vright[jj] -= model.dVleft[j]
            end
            for i in (model.id_params+1):(model.nb_params_maintenance - 1)
                ii = ind_ij(i, model.id_params + 1)
                model.d2Vright[ii] -= model.dVleft[i]
            end
        end
        for i in 1:npm
             for j in 1:i
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                ij = ind_ij(i, j)
                model.d2A[ij] = m.ρ_QR^δ * model.d2A[ij]
            end
        end
        for j in 1:model.id_params
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            jj = ind_ij(model.id_params, j)
            model.d2A[jj] += m.ρ_QR^δ * δ / m.ρ_QR * model.dA[j]
        end
        model.d2A[ind_ij(model.id_params, model.id_params)] += m.ρ_QR^δ * δ / m.ρ_QR * (2 * model.dA[model.id_params] + (δ - 1) /m.ρ_QR * model.A)
        for i in (model.id_params+1):npm
            # //id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ii = ind_ij(i, model.id_params)
            model.d2A[ii] += m.ρ_QR^δ * δ / m.ρ_QR * model.dA[i]
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    model.dVR_prec[k * npm + i] = (1 - m.ρ_ARA) * model.dVR_prec[(k - 1) * npm + i]
                end
            end
            model.dVR_prec[k * npm + model.id_params + 1] -= model.VR_prec[k-1]
        end
        if nk > 0
            for i in npm
                model.dVR_prec[i] = (1 - m.ρ_ARA) * model.dA[i] * model.Δt
                model.dVright[i] = (1 - m.ρ_ARA) * model.dVleft[i]
            end
            model.dVR_prec[model.id_params + 1] -= model.A * model.Δt
            model.dVright[model.id_params + 1] -= model.Vleft
        else
            for i in 1:npm
                model.dVright[i] = (1 - m.ρ_ARA) * model.dVleft[i]
            end
            model.dVright[model.id_params + 1] -= model.Vleft
        end
        for i in 1:npm
            model.dA[i] = m.ρ_QR^ δ *  model.dA[i]
        end
        model.dA[model.id_params] += m.ρ_QR^δ * δ / m.ρ_QR * model.A
    end
    model.Vright = (1 - m.ρ_ARA) * model.Vleft
    if nk > 1
        for k in (nk - 1):-1:1 
            model.VR_prec[k + 1] = (1 - m.ρ_ARA) * model.VR_prec[k]
        end
    end
    if nk > 0 
        model.VR_prec[1] = (1 - m.ρ_ARA) * model.A * model.Δt
    end
    model.A = m.ρ_QR^δ * model.A
end

function update!(m::ARAm, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    inc!(model) #model.k += 1;

    nk = model.k
    if nk > model.mu 
        nk = model.mu
    end
    nk2 = nk
    if nk > m.m-1 
        nk2 = m.m-1
    end

    # //printf("ARAinf k=%d,max_mem=%d, nk=%d\n",model.k,model.max_mem,nk);
    npm = model.nb_params_maintenance
    if hessian
        npm2 = ind_nb(npm)
        if nk > nk2
            for k in (nk - 1):-1:nk2
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2VR_prec[k * npm2 + ij] = model.d2VR_prec[(k - 1) * npm2 + ij]
                    end
                end
            end
        end

        if (model.k > m.m) && (nk2 > 0)
            if nk > nk2
                for i in 1:npm
                     for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2Vright[ij] -= m.ρ * model.d2VR_prec[(nk2-1) * npm2 + ij]
                        model.d2VR_prec[nk2 * npm2 + ij] = (1-m.ρ) * model.d2VR_prec[(nk2-1) * npm2 + ij]
                    end
                end
                for j in 1:model.id_params
                    jj = ind_ij(model.id_params, j)
                    model.d2Vright[jj] -= model.dVR_prec[(nk2-1) * npm + j]
                    model.d2VR_prec[nk2 * npm2 + jj] -= model.dVR_prec[(nk2-1) * npm + j]
                end
                for i in (model.id_params + 1):npm
                    ii = ind_ij(i, model.id_params)
                    model.d2Vright[ii] -= model.dVR_prec[(nk2 - 1) * npm + i]
                    model.d2VR_prec[nk2 * npm2 + ii] -= model.dVR_prec[(nk2-1) * npm + i]
                end
            else
                for i in 1:npm
                     for j in 1:i
                        ij = ind_ij(i, j) 
                        model.d2Vright[ij] -= m.ρ * model.d2VR_prec[(nk2-1) * nbd2 + ij]
                    end
                end
                for j in 1:model.id_params 
                    model.d2Vright[ind_ij(model.id_params, j)] -= model.dVR_prec[(nk2 - 1) * npm + j]
                end
                for i in (model.id_params + 1):npm  
                    model.d2Vright[ind_ij(i, model.id_params)] -= model.dVR_prec[(nk2 - 1) * npm + i]
                end
            end
        end

        if nk2 > 1
            for k in (nk2 - 1):-1:1
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2Vright[ij] -= m.ρ * model.d2VR_prec[(k - 1) * npm2 + ij]
                        model.d2VR_prec[k * npm2 + ij] = (1 - m.ρ) * model.d2VR_prec[(k - 1) * npm2 + ij]
                    end
                end
                for j in 1:model.id_params
                    jj = ind_ij(model.id_params, j)
                    model.d2Vright[jj] -= model.dVR_prec[(k - 1) * npm + j]
                    model.d2VR_prec[k * npm2 + jj] -= model.dVR_prec[(k - 1) * npm + j]
                end
                for i in (model.id_params + 1):npm
                    ii = ind_ij(i, model.id_params)
                    model.d2Vright[ii] -= model.dVR_prec[(k - 1) * npm + i]
                    model.d2VR_prec[k * npm2 + ii] -= model.dVR_prec[(k - 1) * npm + i]
                end
            end
        end
        if nk > 0
            for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    prov = (1 - m.ρ) * model.d2A[ij] * model.Δt
                    model.d2VR_prec[ij] = prov
                    model.d2Vright[ij] += prov
                end
            end
            for j in 1:model.id_params
                jj = ind_ij(model.id_params, j)
                prov = model.dA[j] * model.Δt
                model.d2VR_prec[jj] -= prov
                model.d2Vright[jj] -= prov
            end
            for i in (model.id_params + 1):npm
                ii = ind_ij(i, model.id_params)
                prov = model.dA[i]*model.Δt
                model.d2VR_prec[ii] -= prov
                model.d2Vright[ii] -= prov
            end
        else
           for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2Vright[ij] += (1 - m.ρ) * model.d2A[ij] * model.Δt
                end
            end
            for j in 1:model.id_params
                model.d2Vright[ind_ij(model.id_params, j)] -= model.dA[j] * model.Δt
            end
            for i in (model.id_params + 1):npm
                model.d2Vright[ind_ij(i, model.id_params)] -= model.dA[i] * model.Δt
            end
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    model.dVR_prec[k * npm + i] = model.dVR_prec[(k - 1) * npm + i]
                end
            end
        end
        if model.k > m.m && nk2 > 0 
            if nk > nk2
                for i in 1:npm
                    model.dVright[i] -= m.ρ * model.dVR_prec[(nk2 - 1) * npm + i]
                    model.dVR_prec[nk2 * npm + i] = (1 - m.ρ) * model.dVR_prec[nk2 * npm + i]
                end
                model.dVR_prec[nk2 * npm + model.id_params] -= model.VR_prec[nk2]
                model.dVright[model.id_params] -= model.VR_prec[nk2]
            else
                for i in 1:npm 
                    model.dVright[i] -= m.ρ * model.dVR_prec[(nk2 - 1) * npm + i]
                end
                model.dVright[model.id_params] -= model.VR_prec[nk2 - 1]
            end
        end
        if nk2 > 1
            for k in (nk2 - 1):-1:1
                for i in 1:npm
                    model.dVright[i] -= m.ρ * model.dVR_prec[(k - 1) * npm + i]
                    model.dVR_prec[k * npm + i] = (1 - m.ρ) * model.dVR_prec[(k - 1) * npm + i]
                end
                model.dVR_prec[k * npm + model.id_params] -= model.VR_prec[k-1]
                model.dVright[model.id_params] -= model.VR_prec[k-1]
            end
        end
        if nk > 0
            for i in 1:npm
                prov = (1 - m.ρ) * model.dA[i] * model.Δt
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
            prov = model.A * model.Δt
            model.dVR_prec[model.id_params] -= prov
            model.dVright[model.id_params] -=prov
        else
            for i in 1:npm
                model.dVright[i] += (1 - m.ρ) * model.dA[i] * model.Δt
            end
            model.dVright[model.id_params] -= model.A * model.Δt
        end
    end
    # //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
    if nk > nk2 + 1
        for k=(nk - 1):-1:(nk2 + 1)
            model.VR_prec[k + 1]=model.VR_prec[k]
        end
    end
    if model.k > m.m && nk2 > 0
        if nk > nk2
            model.Vright -= m.ρ * model.VR_prec[nk2]
            model.VR_prec[nk2 + 1] = (1-m.ρ) * model.VR_prec[nk2]
        else 
            model.Vright -= m.ρ * model.VR_prec[nk2]
        end
    end
    if nk2 > 1
        for k in (nk2 - 1):-1:1
            model.Vright -= m.ρ * model.VR_prec[k]
            model.VR_prec[k + 1] = (1 - m.ρ) * model.VR_prec[k]
        end
    end
    prov = (1 - m.ρ) * model.A * model.Δt
    #//printf("Vright=%f, m.ρ=%f, A=%f, Tk=%f, Tk-1=%f\n",model.Vright,m.ρ,model.A,model.time[model.k],model.time[model.k - 1]);
    if nk > 0 
        model.VR_prec[1] = prov
    end
    model.Vright += prov
    #//printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model.Vright,nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
end

function update!(m::GQR_ARAm, model::AbstractModel;gradient::Bool=false,hessian::Bool=false)
    inc!(model) #model.k += 1;
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
    npm = model.nb_params_maintenance
    δ = (m.f(K) - m.f(K-1))
    if hessian
        npm2 = ind_nb(npm)
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2VR_prec[k * npm2 + ij] = model.d2VR_prec[(k - 1) * npm2 + ij]
                    end
                end
            end
        end

        if  (model.k >= m) && (nk2 > 0)
            if nk > nk2
                for i in 1:npm
                    for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2Vright[ij] -= m.ρ_ARA * model.d2VR_prec[(nk2 - 1) * npm2 + ij]
                        model.d2VR_prec[nk2 * npm2 + ij] = (1 - m.ρ_ARA) * model.d2VR_prec[(nk2 - 1) * npm2 + ij]
                    end
                end
                for j in 1:model.id_params
                    jj = ind_ij(model.id_params + 1, j)
                    model.d2Vright[jj] -= model.dVR_prec[(nk2 - 1) * npm + j]
                    model.d2VR_prec[nk2 *  npm2 + jj] -= model.dVR_prec[(nk2 - 1) * npm + j]
                end
                for i in (model.id_params+2):npm
                    ii = ind_ij(i, model.id_params + 1)
                    model.d2Vright[ii] -= model.dVR_prec[(nk2 - 1) * npm + i]
                    model.d2VR_prec[nk2 * npm2 + ii] -= model.dVR_prec[(nk2 - 1) * npm + i]
                end
            else
                for i in 1:npm
                     for j in 1:i
                        ij = ind_ij(i, j)
                        model.d2Vright[ij] -= m.ρ_ARA * model.d2VR_prec[(nk2 - 1) * npm2 + ij]
                    end
                end
                for j in 1:(model.id_params + 1) 
                    model.d2Vright[ind_ij(model.id_params + 1, j)] -= model.dVR_prec[(nk2 - 1) * npm +j]
                end
                for i in (model.id_params + 2):npm 
                    model.d2Vright[ind_ij(i, model.id_params + 1)] -= model.dVR_prec[(nk2 - 1) * npm + i]
                end
            end
        end

        for k in nk2:-1:1
            for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2Vright[ij] -= m.ρ_ARA * model.d2VR_prec[(k - 1) * npm2 + ij]
                    model.d2VR_prec[k * npm2 + ij] = (1 - m.ρ_ARA) * model.d2VR_prec[(k - 1) * npm2 + ij]
                end
            end
            for j in 1:(model.id_params + 1)
                jj = ind_ij(model.id_params + 1, j)
                model.d2Vright[jj] -= model.dVR_prec[(k - 1) * npm + j]
                model.d2VR_prec[k * npm2 + jj] -= model.dVR_prec[(k - 1) * npm + j]
            end
            for i in (model.id_params + 2):npm
                ii = ind_ij(i, model.id_params + 1)
                model.d2Vright[ii] -= model.dVR_prec[(k - 1) * npm + i]
                model.d2VR_prec[k * npm2 + ii] -= model.dVR_prec[(k - 1) * npm + i]
            end
        end
        if nk > 0
            for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    prov = (1 - m.ρ_ARA) * model.d2A[ij] * model.Δt
                    model.d2VR_prec[ij] = prov
                    model.d2Vright[ij] += prov
                end
            end
            for j in 1:(model.id_params + 1)
                prov = model.dA[j] * model.Δt
                jj = ind_ij(model.id_params + 1, j)
                model.d2VR_prec[jj] -= prov
                model.d2Vright[jj] -= prov
            end
            for i in (model.id_params + 2):npm
                prov=model.dA[i] * model.Δt
                ii = ind_ij(i, model.id_params + 1)
                model.d2VR_prec[ii] -= prov
                model.d2Vright[ii] -= prov
            end
        else
           for i in 1:npm
                 for j in 1:i
                    ij = ind_ij(i, j)
                    model.d2Vright[ij] += (1 - m.ρ_ARA) * model.d2A[ij] * model.Δt
                end
            end
            for j in 1:(model.id_params + 1)
                model.d2Vright[ind_ij(model.id_params + 1, j)] -= model.dA[j] * model.Δt
            end
            for i in (model.id_params + 2):npm
                model.d2Vright[ind_ij(i, model.id_params + 1)] -= model.dA[i] * model.Δt
            end
        end

        for i in 1:npm
            for j in 1:i
               ij = ind_ij(i, j)
                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
                model.d2A[ij] =m.ρ_QR^δ * model.d2A[ij]
            end
        end
        for j in 1:(model.id_params + 1)
            #//i(<=model.id_params) and id respectively correspond to the column and line indices of (inferior diagonal part of) the hessian matrice
            model.d2A[ind_ij(model.id_params, j)] += m.ρ_QR^δ * δ / m.ρ_QR * model.dA[j]
        end
        model.d2A[ind_ij(model.id_params, model.id_params)] += m.ρ_QR^δ * δ / m.ρ_QR * (2 * model.dA[model.id_params] + (δ - 1) / m.ρ_QR * model.A)
        for i in (model.id_params + 2):npm
            # //id and i(>=model.id_params) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            model.d2A[ind_ij(i, model.id_params)] += m.ρ_QR^δ * δ / m.ρ_QR * model.dA[i]
        end
    end
    if gradient || hessian
        if nk > 1
            for k in (nk - 1):-1:1
                for i in 1:npm
                    model.dVR_prec[k * npm + i] = model.dVR_prec[(k - 1) * npm + i]
                end
            end
        end
        if  (model.k >= m) && (nk2 > 0)
            if nk > nk2
                for i in 1:npm
                    model.dVright[i] -= m.ρ_ARA * model.dVR_prec[(nk2 - 1) * npm + i]
                    model.dVR_prec[nk2 * npm + i] = (1 - m.ρ_ARA) * model.dVR_prec[(nk2 - 1) * npm + i]
                end
                model.dVR_prec[nk2 * npm + model.id_params + 1] -= model.VR_prec[nk2 - 1]
                model.dVright[model.id_params + 1] -= model.VR_prec[nk2 - 1]
            else
                for i in 1:npm
                    model.dVright[i] -= m.ρ_ARA * model.dVR_prec[(nk2 - 1) * npm + i]
                end
                model.dVright[model.id_params + 1] -= model.VR_prec[nk2 - 1]
            end
        end
        for k in nk2:-1:1
            for i in 1:npm
                model.dVright[i] -= m.ρ_ARA * model.dVR_prec[(k - 1) * npm + i]
                model.dVR_prec[k * npm + i] = (1 - m.ρ_ARA) * model.dVR_prec[(k - 1) * npm + i]
            end
            model.dVR_prec[k * npm + model.id_params + 1] -= model.VR_prec[k - 1]
            model.dVright[model.id_params + 1] -= model.VR_prec[k - 1]
        end
        if nk > 0
            for i in 1:npm
                prov = (1 - m.ρ_ARA) * model.dA[i] * model.Δt
                model.dVR_prec[i] = prov
                model.dVright[i] += prov
            end
            prov = model.A * model.Δt
            model.dVR_prec[model.id_params + 1] -= prov
            model.dVright[model.id_params + 1] -= prov
        else
            for i in 1:npm
                model.dVright[i] += (1 - m.ρ_ARA) * model.dA[i] * model.Δt
            end
            model.dVright[model.id_params + 1] -= model.A * model.Δt
        end
        for i in 1:npm
            model.dA[i] = m.ρ_QR^δ *  model.dA[i]
        end
        model.dA[model.id_params] += m.ρ_QR^δ * δ / m.ρ_QR * model.A
    end
    # //printf("Avant nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
    for k=nk:-1:nk2
        model.VR_prec[k] = model.VR_prec[k-1]
    end
    if (model.k >= m) && (nk2 > 0)
        if nk > nk2
            model.Vright -= m.ρ_ARA * model.VR_prec[nk2 - 1]
            model.VR_prec[nk2] = (1 - m.ρ_ARA) * model.VR_prec[nk2 - 1]
        else 
            model.Vright -= m.ρ_ARA * model.VR_prec[nk2 - 1]
        end
    end
    for k=nk2:-1:1
        model.Vright -= m.ρ_ARA* model.VR_prec[k - 1]
        model.VR_prec[k] = (1 - m.ρ_ARA) * model.VR_prec[k - 1]
    end
    prov = (1 - m.ρ_ARA) * model.A * model.Δt
    #//printf("Vright=%f, m.ρ=%f, A=%f, Tk=%f, Tk-1=%f\n",model.Vright,m.ρ,model.A,model.time[model.k],model.time[model.k - 1]);
    if nk > 0 
        model.VR_prec[1] = prov
    end
    model.Vright += prov
    #//printf("Apres Vright=%f, nk=%i, nk2=%i, VRprec[0]=%f, VRprec[1]=%f, VRprec[2]=%f\n",model.Vright,nk,nk2,model.VR_prec[0],model.VR_prec[1],model.VR_prec[2]); 
    model.A = m.ρ_QR^δ * model.A
end