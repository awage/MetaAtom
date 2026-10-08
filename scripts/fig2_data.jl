# Data for Fig. 2: number of attractors, basin entropy and boundary basin
# entropy on a (δ, σ) grid. Run with several threads: julia -t auto fig2_data.jl
#
# Each basin is cached by get_basins, so an interrupted run resumes where it
# stopped. Parameter points that fail (an exception, or initial conditions
# labelled -1 because the mapper did not converge) get NaN and are listed in
# data/fig2_failed.tsv. To redo them: tighten the tolerances in get_mapper /
# get_smap, set redo_failed = true and run again. Only the listed files are
# recomputed; everything else is loaded from the cache.
using DrWatson
@quickactivate
using Attractors
using ProgressMeter

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

redo_failed = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
res = 500
yg = range(-5,5, length=res); grid = (yg, yg)
Np = 200
δrange = range(-30, 30, length = Np)
σrange = range(0.1, 0.5, length = Np)
failed_file = datadir("fig2_failed.tsv")

to_redo = Set{String}()
if redo_failed && isfile(failed_file)
    to_redo = Set(first(split(l, '\t')) for l in readlines(failed_file)[2:end])
    println("Recomputing $(length(to_redo)) failed basins")
end

Sb = fill(NaN, Np, Np)
Sbb = fill(NaN, Np, Np)
Na = fill(NaN, Np, Np)
failed = Tuple{String, Float64, Float64, String}[]
lk = ReentrantLock()
cases = vec(CartesianIndices((Np, Np)))
prog = Progress(length(cases))
Threads.@threads :dynamic for c in cases
    j, k = Tuple(c)
    dp = deepcopy(dps)
    dp.δ = δrange[j]
    dp.σ = σrange[k]
    file = basins_file(dp, grid)
    reason = ""
    try
        data = get_basins(dp, grid; force = file in to_redo, show_progress = false)
        @unpack bas = data
        nlost = count(==(-1), bas)
        if nlost > 0
            reason = "$nlost initial conditions labelled -1"
        else
            Sb[j,k], Sbb[j,k] = basin_entropy(bas, 10)
            Na[j,k] = length(unique(bas))
        end
    catch e
        e isa InterruptException && rethrow()
        reason = replace(sprint(showerror, e), r"\s+" => " ")
    end
    if !isempty(reason)
        lock(lk) do
            push!(failed, (file, dp.δ, dp.σ, reason))
        end
    end
    next!(prog)
end
finish!(prog)

open(failed_file, "w") do io
    println(io, "file\tdelta\tsigma\treason")
    for (file, δ, σ, reason) in sort(failed)
        println(io, file, '\t', δ, '\t', σ, '\t', reason)
    end
end
println("$(length(failed)) of $(length(cases)) parameter points failed, listed in $failed_file")

wsave(datadir("fig2_data.jld2"), @strdict Sb Sbb Na δrange σrange dps grid)
