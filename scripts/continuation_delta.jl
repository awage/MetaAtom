using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

function compute_delta_sweep(params::Dict)
    @unpack  δrange, Np, Nsamples, dps, grid= params


    # compute globl continuation of attractors
    sampler, _ = statespace_sampler(grid) 
    mapper = get_mapper(dps)

    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pidx = :δ

    fractions_cont, attractors_cont = global_continuation(
        ascm, δrange, pidx, sampler; samples_per_parameter = Nsamples
    )

    branches = get_branches(δrange, pidx, dps, attractors_cont)

    return @strdict(fractions_cont,  attractors_cont, δrange, dps, branches)
end

force = true
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
res = 100
yg = range(-5,5, length=res); grid = (yg, yg) 
Np = 15; Nsamples = 10
δrange = range(-30, 30, length = Np)


# Compute the fractions and attractor branches
params = @strdict Np Nsamples δrange dps grid
dat, _ = produce_or_load(compute_delta_sweep, params, datadir(); prefix = "delta_sweep", force)

# Compute the basin entropy for the same range of parameters
dat_ent = get_entropy_δ_sweep(δrange, dps, grid; force)



# Get variables from dictionary
@unpack fractions_cont, attractors_cont, branches = dat

# First paint the bands for the volume of each attractors
colors = colors_from_keys(keys(fractions_cont))
fig = plot_basins_curves(fractions_cont, δrange; colors)

# Add the bifurcation diagram.

# ss = continuation_series(fraction_cont)
ax = Axis(fig[2,1], ylabel = "xn", yticklabelsize = 10, xticklabelsize = 20, ylabelsize = 20)
for k in keys(branches)
    P = StateSpaceSet(branches[k])
    scatter!(ax, P[:,1],P[:,2], markersize = 2.1, color = colors[k], rasterize = false)
end
xlims!(ax,δrange[1],δrange[end])

# Add the basin entropy
@unpack Sb, Sbb = dat_ent
ax = Axis(fig[3,1], ylabel = "Sb", yticklabelsize = 10, xticklabelsize = 20, ylabelsize = 20, xlabel = L"\delta")
scatter!(ax, δrange, Sb; markersize = 5)
xlims!(ax,δrange[1],δrange[end])

save("test.pdf",fig)
