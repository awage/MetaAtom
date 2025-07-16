using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

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

    branches = get_branches(prange, pidx, dps, attractors_cont)

    return @strdict(fractions_cont,  attractors_cont, prange, dps, branches)
end


Np = 250; Nsamples = 5000
params = @strdict Np Nsamples

dat, _ = produce_or_load(compute_delta_sweep, params; prefix = "delta_sweep", force = false)

@unpack fractions_cont, attractors_cont, prange, branches = dat



colors = colors_from_keys(keys(fractions_cont))
fig = plot_basins_curves(fractions_cont, prange; colors)

# fig = plot_attractors_curves(fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange)
#
# P = get_bif_points(branches)
# fig = Figure(size = (1524, 568))

ax = Axis(fig[2,1], ylabel = "xn", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20, xlabel = L"\delta")
for (j,p) in enumerate(branches)
    P = StateSpaceSet(p[2])
    scatter!(ax, P[:,1],P[:,2], markersize = 2.7, color = colors[j], rasterize = true)
end
