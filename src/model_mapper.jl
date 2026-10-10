using Attractors
# using OrdinaryDiffEq:Vern9
using OrdinaryDiffEqVerner
using LinearAlgebra
using SHA


mutable struct MetaAtomParameters{N}
    ω::N
    σ::N
    β::N
    η::N
    μ::N
    δ::N
end

# Cache file names: SHA-1 of the exact parameter values. Unlike savename, no
# rounding (δ = -25.99 and -25.98 get different files) and nothing is dropped
# (grids and ranges are part of the key). Stable across Julia sessions.
canonical(p::MetaAtomParameters) = (p.ω, p.σ, p.β, p.η, p.μ, p.δ)
canonical(r::AbstractRange) = (first(r), last(r), length(r))
canonical(t::Tuple) = map(canonical, t)
canonical(x) = x
hash_name(args...) = bytes2hex(sha1(repr(map(canonical, args))))[1:16]

    
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
    yg =  collect(range(-30, 30; length = 1501))
    grid = (yg, yg)
    mapper = AttractorsViaRecurrences(smap, grid; 
                    consecutive_basin_steps = 10, 
                    consecutive_recurrences = 2000,
                    attractor_locate_steps = 1000,
                    maximum_iterations = 1e5)
    return mapper
end

function get_smap(dps::MetaAtomParameters)
    (;ω, σ, β, η, μ, δ) = dps
    diffeq = (alg = Vern9(), reltol = 1e-8, maxiters = 1e10)
    ds = CoupledODEs(oscillator_model!, rand(2), dps; diffeq)
    smap = StroboscopicMap(ds, 2*pi/ω) # Stroboscopic map definition
    return smap
end

function get_trajectory_map(u0, T, dps::MetaAtomParameters)
    smap = get_smap(dps)
    y,t = trajectory(smap, T, u0) 
    return  y,t 
end

function get_trajectory_ode(u0, T, dps::MetaAtomParameters)
    (;ω, σ, β, η, μ, δ) = dps
    diffeq = (alg = Vern9(), reltol = 1e-8, maxiters = 1e6)
    ds = CoupledODEs(oscillator_model!, rand(2), dps; diffeq)
    y,t = trajectory(ds, T, u0; Δt = 0.01) 
    return  y,t 
end
