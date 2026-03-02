using DrWatson
@quickactivate 
using CairoMakie
using ProgressMeter
include(srcdir("model_mapper.jl"))


function compute_basins(mapper, xg, yg)
    bsn = @showprogress [ mapper([x; y]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    return bsn,  att
end

res = 300

f = Figure(size = (400,400))
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = -25.9; β = 0.4;
dps = model_parameters(ω, σ, β, η, μ, δ)
yg = range(-80,40, length = res) 
xg = range(-15,25, length = res) 
mapper = get_mapper(dps)
bsn, att = compute_basins(mapper, xg, yg)
labs_args = (ylabel = L"i_0", xlabel = L"q_0", yticklabelsize = 20,  ylabelsize = 25, xticklabelsize = 20,  xlabelsize = 25)
ax = Axis(f[1,1]; labs_args...) #, yscale = log10);
heatmap!(ax, yg, xg, bsn'; rasterize = true, colormap = ColorScheme([RGB(0.1,0.1,0.1), RGB(1,0.1,0.1),  RGB(0.24,0.24,1)]))
save("basins_delta_-25.9.pdf",f)

f = Figure(size = (400,400))
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = -25.0; β = 0.4; 
dps = model_parameters(ω, σ, β, η, μ, δ)
yg = range(-80,40, length = res) 
xg = range(-15,25, length = res) 
mapper = get_mapper(dps)
bsn, att = compute_basins(mapper, xg, yg)
labs_args = (ylabel = L"i_0", xlabel = L"q_0", yticklabelsize = 20,  ylabelsize = 25, xticklabelsize = 20,  xlabelsize = 25)
ax = Axis(f[1,1]; labs_args...) #, yscale = log10);
heatmap!(ax, yg, xg, bsn'; rasterize = true, colormap = ColorScheme([RGB(0.1,0.1,0.1), RGB(1,0.1,0.1),  RGB(0.24,0.24,1)]))
save("basins_delta_-25.0.pdf",f)

