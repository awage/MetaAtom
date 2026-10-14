# Transient lengths on the (δ, σ) grid of Fig. 2. Run with several threads:
#   julia -t auto fig2_transients.jl
#
# For every parameter point the attractors are loaded from the basins cached by
# fig2_data.jl: the fine grid file if it exists, else the coarse grid file (the
# fine pass only ran where the coarse pass found two attractors or more). A
# BasinMapProximity with these attractors maps the res_ic² initial conditions,
# and convergence_time gives the number of map periods (T = 2π/ω) each one takes
# to come ε-close to its attractor.
#
# Fallbacks:
# - ε: the automatic value of BasinMapProximity, capped at ε_max so that the
#   transients of different points are comparable. If it cannot be computed
#   (a single attractor made of one point), ε_max.
# - ICs that do not converge within maxsteps map periods are mapped again with
#   the recurrences mapper of get_mapper, seeded with the same attractors; its
#   step count is used. If this fails too, the IC is counted as lost.
# - Points with no cached basin (failed in fig2_data.jl) get NaN.
#
# Every point is cached in data/transients, so an interrupted run resumes where
# it stopped. Points that fail are listed in data/fig2_transients_failed.tsv.
using DrWatson
@quickactivate
using Attractors
using ProgressMeter
using Statistics

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
grid_of(res) = (range(-5, 5, length = res), range(-5, 5, length = res))
grid_coarse = grid_of(50)
grid_fine = grid_of(500)
Np = 200
δrange = range(-30, 30, length = Np)
σrange = range(0.1, 0.5, length = Np)
all_cases = vec(CartesianIndices((Np, Np)))

res_ic = 10
grid_ic = grid_of(res_ic)
ε_max = 0.05
maxsteps = 10^6
force = false
failed_file = datadir("fig2_transients_failed.tsv")

# Attractors of dp from the cached basins, and the grid they come from
function cached_attractors(dp)
    for (g, src) in ((grid_fine, "fine"), (grid_coarse, "coarse"))
        file = basins_file(dp, g)
        isfile(file) || continue
        att = try
            wload(file)["att"]
        catch e
            e isa InterruptException && rethrow()
            continue
        end
        att = filter(a -> length(a[2]) > 0, att)
        isempty(att) || return att, src
    end
    return nothing, "none"
end

function proximity_mapper(dp, att)
    ds = get_smap(dp)
    ε = try
        BasinMapProximity(ds, att).ε
    catch e
        e isa InterruptException && rethrow()
        ε_max
    end
    ε = ε > 0 ? min(ε, ε_max) : ε_max
    return BasinMapProximity(ds, att; ε, consecutive_lost_steps = maxsteps), ε
end

function _transients(d)
    @unpack dp, grid_ic = d
    att, src = cached_attractors(dp)
    isnothing(att) && error("no cached basin")
    bmap, ε = proximity_mapper(dp, att)
    ics = vec([[x, y] for x in grid_ic[1], y in grid_ic[2]])
    labels = zeros(Int, length(ics))
    T = fill(Inf, length(ics))
    for (i, u0) in enumerate(ics)
        labels[i] = bmap(u0)
        T[i] = convergence_time(bmap)
    end
    # Fallback for the ICs that did not converge
    redo = findall(isinf, T)
    nfallback = length(redo)
    if nfallback > 0
        try
            mapper = get_mapper(dp; att)
            for i in redo
                labels[i] = mapper(ics[i]; show_progress = false)
                labels[i] > 0 && (T[i] = mapper.bsn_nfo.safety_counter)
            end
        catch e
            e isa InterruptException && rethrow()
            @warn "recurrences fallback failed at δ = $(dp.δ), σ = $(dp.σ)" exception = e
        end
    end
    return @strdict T labels ε src nfallback
end

function get_transients(dp; force = false)
    data, _ = produce_or_load(_transients, @strdict(dp, grid_ic), datadir("transients");
        filename = _ -> hash_name(dp, grid_ic, ε_max, maxsteps),
        prefix = "model_transients", storepatch = false, suffix = "jld2", force)
    return data
end

Tmax = fill(NaN, Np, Np)
Tmean = fill(NaN, Np, Np)
Tmedian = fill(NaN, Np, Np)
Nlost = fill(NaN, Np, Np)
Nfallback = fill(NaN, Np, Np)
εs = fill(NaN, Np, Np)
fine_source = falses(Np, Np)
failed = Tuple{Float64, Float64, String}[]
lk = ReentrantLock()
prog = Progress(length(all_cases); desc = "Transients: ")
# :greedy hands out one case at a time: the cost per case is very uneven
Threads.@threads :greedy for c in all_cases
    j, k = Tuple(c)
    dp = deepcopy(dps)
    dp.δ = δrange[j]
    dp.σ = σrange[k]
    try
        @unpack T, ε, src, nfallback = get_transients(dp; force)
        conv = filter(isfinite, T)
        Nlost[j,k] = length(T) - length(conv)
        Nfallback[j,k] = nfallback
        εs[j,k] = ε
        fine_source[j,k] = src == "fine"
        if !isempty(conv)
            Tmax[j,k] = maximum(conv)
            Tmean[j,k] = mean(conv)
            Tmedian[j,k] = median(conv)
        end
    catch e
        e isa InterruptException && rethrow()
        lock(() -> push!(failed, (dp.δ, dp.σ, replace(sprint(showerror, e), r"\s+" => " "))), lk)
    end
    next!(prog)
end
finish!(prog)

open(failed_file, "w") do io
    println(io, "delta\tsigma\treason")
    for (δ, σ, reason) in sort(failed)
        println(io, δ, '\t', σ, '\t', reason)
    end
end
println("$(length(failed)) of $(length(all_cases)) points failed, listed in $failed_file")
println("$(count(>(0), filter(!isnan, Nfallback))) points used the recurrences fallback, ",
    "$(count(>(0), filter(!isnan, Nlost))) have lost ICs")

wsave(datadir("fig2_transients.jld2"), @strdict Tmax Tmean Tmedian Nlost Nfallback εs fine_source δrange σrange dps grid_ic ε_max maxsteps)
