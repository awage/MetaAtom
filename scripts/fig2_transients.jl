using DrWatson
@quickactivate
using Attractors
using ProgressMeter

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

redo_failed = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
res_coarse = 50
res_fine = 500
grid_of(res) = (range(-5, 5, length = res), range(-5, 5, length = res))
grid_coarse = grid_of(res_coarse)
grid = grid_of(res_fine)
Np = 200
δrange = range(-30, 30, length = Np)
σrange = range(0.1, 0.5, length = Np)
force = false
all_cases = vec(CartesianIndices((Np, Np)))


transient_time(m) = m.bsn_nfo.safety_counter * m.bsn_nfo.Δt

res = 10
grid_transient = grid_of(res)
transients = zeros(Np,Np)

for c in all_cases
    j, k = Tuple(c)
    dp = deepcopy(dps)
    dp.δ = δrange[j]
    dp.σ = σrange[k]
    data = get_basins(dp, grid_coarse; force = false, show_progress = false)
    @unpack att = data
    mapper = get_mapper(dp)
    # Seed the mapper with the known attractors (a = label => StateSpaceSet)
    for (_, A) in att; mapper(A[1]); end
    t_tmp = Float64[]
    for x in grid_transient[1], y in grid_transient[2]
        mapper([x, y])
        push!(t_tmp, transient_time(mapper))
    end
    transients[j,k] = maximum(t_tmp)
end

wsave(datadir("fig2_transients.jld2"), @strdict transients grid_transients)


