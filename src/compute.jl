mutable struct Compute
    S1::Float64
    S2::Float64
    S0::Float64
    S3::Float64
    S4::Float64

    dS1::Vector{Float64}
    dS2::Vector{Float64}
    dS3::Vector{Float64}
    dS4::Vector{Float64} 

    d2S1::Vector{Float64}
    d2S2::Vector{Float64}
    d2S3::Vector{Float64}
end

function init!(c::Compute, m::AbstractModel; with_deriv::Bool = false)
    c.S0, c.S1, c.S2, c.S3, c.S4 = zeros(5)
    if with_deriv
        nbd = m.nb_params_maintenance + m.nb_params_family - 1
        nb2d = nbd * (nbd + 1) ÷ 2
        c.dS1, c.dS2, c.dS3=zeros(nbd), zeros(nbd), zeros(m.nb_params_maintenance)
        if m.nb_params_cov > 0
            c.dS4 = zeros(m.nb_params_cov)
        end
        c.d2S1, c.d2S2 = zeros(nb2d), zeros(nb2d)  #inferior diagonal part of the hessian matrice by lines
        c.d2S3 = zeros(m.nb_params_maintenance * (nb_params_maintenance + 1) ÷ 2)
    end
end