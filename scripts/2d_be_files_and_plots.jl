using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using ProgressMeter

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

force = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
res = 500
yg = range(-5,5, length=res); grid = (yg, yg) 
Np = 100; 
δrange = range(-30, 30, length = Np)
σrange = range(0.1, 0.5, length = Np)

Sb = zeros(Np, Np)
Sbb = zeros(Np, Np)
Na = zeros(Np, Np)
@showprogress for j in eachindex(δrange) 
    for k in eachindex(σrange)
        # @show δrange[j], σrange[k]
        dp = deepcopy(dps) 
        dp.δ = δrange[j]
        dp.σ = σrange[k]
        try
            data = get_basins_tmp(dp, grid; force, show_progress = true)
            @unpack bas = data
            Sb[j,k], Sbb[j,k] = basin_entropy(bas, 10) 
            Na[j,k] = length(unique(bas))
        catch 
            Sb[j,k] = Sbb[j,k] = 0.0
        end
    end
end

using JLD2 
@save "basin_entropy_model.jld2" Sb Sbb Na δrange σrange dps
fig = Figure(size = (600,600))
ax = Axis(fig[1,1], ylabel = L"\delta", xlabel = L"\sigma",  yticklabelsize = 10, xticklabelsvisible = false, ylabelsize = 15, xlabelsize = 15)
heatmap!(ax, δrange, σrange, Sb; rasterize = true)
save("2d_basin_entropy.pdf", fig)

fig = Figure(size = (600,600))
ax = Axis(fig[1,1], ylabel = L"\delta", xlabel = L"\sigma",  yticklabelsize = 10, xticklabelsvisible = false, ylabelsize = 15, xlabelsize = 15)
heatmap!(ax, δrange, σrange, Na; rasterize = true)
save("2d_basins_Na.pdf", fig)

fig = Figure(size = (600,600))
ax = Axis(fig[1,1], ylabel = L"\delta", xlabel = L"\sigma",  yticklabelsize = 10, xticklabelsvisible = false, ylabelsize = 15, xlabelsize = 15)
heatmap!(ax, δrange, σrange, Sbb; rasterize = true)
save("2d_basins_Sbb.pdf", fig)
