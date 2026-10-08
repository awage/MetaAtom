using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using ProgressMeter

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4

using JLD2 
@load "basin_entropy_model_Np=200_verif.jld2"
# @load "basin_entropy_model.jld2"
using ColorSchemes
# for sch in keys(colorschemes)

sch = :jet1 
args =  (yticklabelsize = 25, xticklabelsvisible = true, ylabelsize = 35, xticklabelsize = 25, xlabelsize = 35)

    fig = Figure(size = (600,600))
    ax = Axis(fig[1,1];  xlabel = L"\delta", ylabel = L"\sigma", args... )
    heatmap!(ax, δrange, σrange, Sb; rasterize = true, colormap = sch)
    Colorbar(fig[1, 2], limits = (0, maximum(Sb)), colormap = sch, ticklabelsize = 25)
    save(string("fig2b",  ".png"), fig)

    fig = Figure(size = (600,600))
    ax = Axis(fig[1,1];  xlabel = L"\delta", ylabel = L"\sigma",  args...)
    heatmap!(ax, δrange, σrange, Sbb; rasterize = true, colormap = sch)
    Colorbar(fig[1, 2],  limits = (0, maximum(filter(!isnan,Sbb))), colormap = sch, ticklabelsize = 25)
    save(string("fig2c",  ".png"), fig)

    fig = Figure(size = (600,600))
    ax = Axis(fig[1,1];  xlabel = L"\delta", ylabel = L"\sigma",  args...)
    heatmap!(ax, δrange, σrange, Na; rasterize = true, colormap = sch)
    Colorbar(fig[1, 2],  limits = (1, maximum(Na)), colormap = sch, ticklabelsize = 25)
    save(string("fig2a",  ".png"), fig)
# end
