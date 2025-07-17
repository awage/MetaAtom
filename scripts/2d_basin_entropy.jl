using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

force = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
res = 100
yg = range(-5,5, length=res); grid = (yg, yg) 
Np = 5; 
δrange = range(-30, 30, length = Np)
σrange = range(0.1, 0.5, length = Np)

Sb = zeros(Np, Np)
Sbb = zeros(Np, Np)
for j in eachindex(δrange) 
    @Threads.threads for k in eachindex(σrange)
        @show δrange[j], σrange[k]
        dp = deepcopy(dps) 
        dp.δ = δrange[j]
        dp.σ = σrange[k]
        data = get_basins(dp, grid; force, show_progress = false)
        @unpack bas = data
        Sb[j,k], Sbb[j,k] = basin_entropy(bas, 10) 
    end
end
