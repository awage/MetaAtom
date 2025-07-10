using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("model_mapper.jl"))


function compute_delta_sweep(params::Dict)
    @unpack  Np, Nsamples = params

    σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
    prange = range(-30, 30, length = Np)
    dps = model_parameters(ω, σ, β, η, μ, δ)

    # compute globl continuation of attractors
    yg = range(-5,5, length=10)
    grid = (yg, yg) 
    sampler, _ = statespace_sampler(grid) 
    mapper = get_mapper(dps)

    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pidx = :δ

    fractions_cont, attractors_cont = global_continuation(
        ascm, prange, pidx, sampler; samples_per_parameter = Nsamples
    )

    return @strdict(fractions_cont,  attractors_cont, prange, dps)
end


Np = 500; Nsamples = 1000
params = @strdict Np Nsamples

dat, _ = produce_or_load(compute_delta_sweep, params; prefix = "delta_sweep", force = false)

@unpack fractions_cont, attractors_cont, prange = dat

fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange,
)

