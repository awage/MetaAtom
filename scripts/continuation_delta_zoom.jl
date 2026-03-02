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
colors = colors_from_keys(unique_keys(fractions_cont))

lab_args = (yticklabelsize = 10, xticklabelsize = 10,   ylabelsize = 17, xlabelsize = 17, xminorticks = IntervalsBetween(5), xminorticksvisible = true, xminorgridvisible = true)

fig = Figure(size = (400, 400))
ax = Axis(fig[1,1]; ylabel = "Fractions", xticklabelsvisible = false,lab_args...)
Attractors.plot_basins_curves!(ax, fractions_cont, δrange; colors) 

# Custom modification of fig

# Add the bifurcation diagram.

# ss = continuation_series(fraction_cont)
ax = Axis(fig[2,1]; ylabel = L"q",xticklabelsvisible = false, lab_args...)

for k in keys(branches)
    P = StateSpaceSet(branches[k])
    scatter!(ax, P[:,1],P[:,3], markersize = 1.0, color = colors[k], rasterize = false)
end
xlims!(ax,δrange[1],δrange[end])

periods = [ Vector{Vector{Float64}}() for k in 1:length(branches)]
for (j,aa) in enumerate(attractors_cont)
    dp = deepcopy(dps) 
    dp.δ = δrange[j]
    smap = get_smap(dp)
    # p = zeros(Int,length(aa))
    for (k,att) in enumerate(aa)
        l = lyapunov(smap, 1000, u0 = att[2][1]; Ttr = 100)
        # @show length(att[2]), sign(l), att[1] 
       if l < 0 
           push!(periods[att[1]], [δrange[j]; length(att[2])])
       else
           push!(periods[att[1]], [δrange[j]; 32])
       end
    end
end

yticks = ([1, 2, 4, 8, 16, 32], ["1", "2", "4", "8", "16", "Chaos"])
ax = Axis(fig[3,1]; ylabel = "periods",  yticks, xticklabelsvisible = false,yscale = log2, lab_args...)
for k in 1:length(periods)
    P = StateSpaceSet(periods[k])
    scatter!(ax, P[:,1],P[:,2], markersize = 3.7, color = colors[k], rasterize = false)
    lines!(ax, P[:,1],P[:,2]; linewidth = 1, color = colors[k], rasterize = false)
end
xlims!(ax,δrange[1],δrange[end])

# Add the basin entropy
@unpack Sb, Sbb = dat_ent
ax = Axis(fig[4,1]; ylabel = L"S_b",  xlabel = L"\delta", xticklabelsvisible = true, lab_args...)
lines!(ax, δrange, Sb; linewidth = 0.5)
xlims!(ax,δrange[1],δrange[end])
ylims!(ax, 0.0, 0.3)
lines!(ax, [-25; -25], [0.; 0.3]; linewidth = 1.5, linestyle = :dash, color = :black)
lines!(ax, [-22.14 ; -22.14], [0.; 0.3]; linewidth = 1.5, linestyle = :dash, color = :black)
lines!(ax, [-27.47; -27.47], [0.;  0.3]; linewidth = 1.5, linestyle = :dash, color = :black)
lines!(ax, [-25.9; -25.9], [0.;  0.3]; linewidth = 1.5, linestyle = :dashdot, color = :black)
lines!(ax, [-28.84;  -28.84], [0.; 0.3]; linewidth = 1.5, linestyle = :dash, color = :black)

save("continuation_zoom.png",fig)
save("continuation_zoom.pdf",fig)
