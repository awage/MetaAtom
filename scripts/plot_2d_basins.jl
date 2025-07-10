using DrWatson
@quickactivate 
using CairoMakie
using ProgressMeter
include(srcdir("model_mapper.jl"))


function compute_basins(mapper)
    xg = yg = range(-5,5, length = 100) 
    bsn = @showprogress [ mapper([x; y]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    return bsn,  att
end

σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)


# compute basins
mapper = get_mapper(dps)
bsn, att = compute_basins(mapper)

xg = yg = range(-5,5, length = 100) 

f = Figure(size = (800,800))
ax = Axis(f[1,1], xlabel = L"x_1", ylabel = L"x_2") #, yscale = log10);

heatmap!(ax, xg, yg, bsn; rasterize = true)

save("test.pdf",f)

