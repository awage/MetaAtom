using Attractors
# using OrdinaryDiffEq:Vern9
using OrdinaryDiffEqVerner
using LinearAlgebra 


mutable struct MetaAtomParameters{N}
    ω::N
    σ::N
    β::N
    η::N
    μ::N
    δ::N
end

    
function model_parameters(ω, σ, β, η, μ, δ)
    return MetaAtomParameters(ω, σ, β, η, μ, δ)
end

function oscillator_model!(du, u, p, t)
    (;ω, σ, β, η, μ, δ) = p
    u1 = u[1]
    u2 = u[2]

    du[1] = u2
    du[2] = -σ*u2 -u1 +β*u1^2 -η*u1^3 +  μ*cos(ω*t) + δ
    return nothing
end

    # du1 = u[2]
    # du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)

function get_mapper(dps::MetaAtomParameters)
    smap = get_smap(dps)
    yg =  collect(range(-20, 20; length = 50001))
    grid = (yg, yg)
    mapper = AttractorsViaRecurrences(smap, grid; 
                    consecutive_basin_steps = 100, 
                    consecutive_recurrences = 1000,
                    attractor_locate_steps = 1000)
    return mapper
end

function get_smap(dps::MetaAtomParameters)
    (;ω, σ, β, η, μ, δ) = dps
    diffeq = (alg = Vern9(), reltol = 1e-8, maxiters = 1e6)
    ds = CoupledODEs(oscillator_model!, rand(2), dps; diffeq)
    smap = StroboscopicMap(ds, 2*pi/ω) # Stroboscopic map definition
    return smap
end

function get_trajectory(u0, T, di::Dict)
    @unpack  N, c, k1, k3, F, kc, ω = di
    W = coupling_matrix(N)
    diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    ds = CoupledODEs(coupled_duffings!, rand(N*2), (N, W, c, k1, k3, F, kc, ω); diffeq)
    smap = StroboscopicMap(ds, 2*pi/ω) # Stroboscopic map definition
    y,t = trajectory(smap, T, u0) 
    return  y,t 
end
