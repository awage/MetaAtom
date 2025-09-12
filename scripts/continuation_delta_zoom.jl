using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using ChaosTools

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

force = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
res = 400
yg = range(-5,5, length=res); grid = (yg, yg) 
Np = 800; Nsamples = 3000
δrange = range(-30, -20, length = Np)


# Compute the fractions and attractor branches
params = @strdict Np Nsamples δrange dps grid
dat, _ = produce_or_load(compute_delta_sweep, params, datadir(); prefix = "delta_sweep_zoom", force)

# Compute the basin entropy for the same range of parameters
dat_ent = get_entropy_δ_sweep(δrange, dps, grid; force = false)

# Get variables from dictionary
@unpack fractions_cont, attractors_cont, branches = dat

# First paint the bands for the volume of each attractors
colors = colors_from_keys(keys(fractions_cont))
fig = plot_basins_curves(fractions_cont, δrange; colors)

# Custom modification of fig
fig.current_axis.x.xticklabelsvisible = false
fig.current_axis.x.xlabel = ""
fig.current_axis.x.ylabel = "Fractions"
fig.current_axis.x.yticklabelsize = 10
fig.current_axis.x.ylabelsize = 15

# Add the bifurcation diagram.

# ss = continuation_series(fraction_cont)
ax = Axis(fig[2,1], ylabel = "xn", yticklabelsize = 10, xticklabelsvisible = false, ylabelsize = 15)
for k in keys(branches)
    P = StateSpaceSet(branches[k])
    scatter!(ax, P[:,1],P[:,3], markersize = 1.5, color = colors[k], rasterize = false)
end
xlims!(ax,δrange[1],δrange[end])

periods = [ Vector{Vector{Float64}}() for k in 1:length(branches)]
for (j,aa) in enumerate(attractors_cont)
    dp = deepcopy(dps) 
    dp.δ = δrange[j]
    smap = get_smap(dp)
    # p = zeros(Int,length(aa))
    for (k,att) in enumerate(aa)
        l = lyapunov(smap, 1000, u0 = att[2][1])
        # @show length(att[2]), sign(l), att[1] 
       if l < 0 
           push!(periods[att[1]], [δrange[j]; length(att[2])])
       else
           push!(periods[att[1]], [δrange[j]; 32])
       end
    end
end

yticks = ([1, 2, 4, 8, 16, 32], ["1", "2", "4", "8", "16", "Chaos"])
ax = Axis(fig[3,1]; ylabel = "periods", yticklabelsize = 10, xticklabelsvisible = false, ylabelsize = 15, yticks, yscale = log2)
for k in 1:length(periods)
    P = StateSpaceSet(periods[k])
    scatter!(ax, P[:,1],P[:,2], markersize = 3.7, color = colors[k], rasterize = false)
    lines!(ax, P[:,1],P[:,2]; linewidth = 1, color = colors[k], rasterize = false)
end
xlims!(ax,δrange[1],δrange[end])

# Add the basin entropy
@unpack Sb, Sbb = dat_ent
ax = Axis(fig[4,1], ylabel = "Sb", yticklabelsize = 10, xticklabelsize = 10, ylabelsize = 15, xlabel = L"\delta", xlabelsize = 20)
scatter!(ax, δrange, Sb; markersize = 5)
xlims!(ax,δrange[1],δrange[end])

save("test.png",fig)
