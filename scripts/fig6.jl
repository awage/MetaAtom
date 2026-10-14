using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("model_mapper.jl"))

using JLD2
@load "../data/fig2_transients.jld2"
using ColorSchemes

# Transients in map periods (T = 2π/ω), from fig2_transients.jl. They span
# several decades, hence the log color scale.
sch = :jet1
args =  (yticklabelsize = 25, xticklabelsvisible = true, ylabelsize = 35, xticklabelsize = 25, xlabelsize = 35)

    fig = Figure(size = (600,600))
    ax = Axis(fig[1,1];  xlabel = L"\delta", ylabel = L"\sigma", args... )
    hm = heatmap!(ax, δrange, σrange, Tmean; rasterize = true, colormap = sch, colorscale = log10)
    Colorbar(fig[1, 2], hm, ticklabelsize = 25)
    save(string("fig6a",  ".png"), fig)

    fig = Figure(size = (600,600))
    ax = Axis(fig[1,1];  xlabel = L"\delta", ylabel = L"\sigma",  args...)
    hm = heatmap!(ax, δrange, σrange, Tmax; rasterize = true, colormap = sch, colorscale = log10)
    Colorbar(fig[1, 2], hm, ticklabelsize = 25)
    save(string("fig6b",  ".png"), fig)
