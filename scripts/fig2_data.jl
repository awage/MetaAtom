# Data for Fig. 2: number of attractors, basin entropy and boundary basin
# entropy on a (δ, σ) grid. Run with several threads: julia -t auto fig2_data.jl
#
# Two passes. The first computes every basin on a coarse res_coarse² grid of
# initial conditions, only to count the attractors. The second recomputes on the
# fine res_fine² grid only the points where the first pass found at least two
# attractors (or failed). Points with a single attractor get Na = 1 and
# Sb = Sbb = 0 without a fine basin; these are the expensive ones when the
# transients are long.
#
# Each basin is cached by get_basins (the grid is part of the file name), so an
# interrupted run resumes where it stopped. Parameter points that fail in the
# fine pass (an exception, or initial conditions labelled -1 because the mapper
# did not converge) get NaN and are listed in data/fig2_failed.tsv. To redo
# them: tighten the tolerances in get_mapper / get_smap, set redo_failed = true
# and run again. Only the listed files are recomputed; everything else is loaded
# from the cache.
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
failed_file = datadir("fig2_failed.tsv")

to_redo = Set{String}()
if redo_failed && isfile(failed_file)
    to_redo = Set(first(split(l, '\t')) for l in readlines(failed_file)[2:end])
    println("Recomputing $(length(to_redo)) failed basins")
end

# Basin of the point (j, k) on grid g, and its Sb, Sbb and Na. A failure is
# returned as a non empty reason, with NaN values.
function basin_stats(j, k, g)
    dp = deepcopy(dps)
    dp.δ = δrange[j]
    dp.σ = σrange[k]
    file = basins_file(dp, g)
    try
        data = get_basins(dp, g; force = file in to_redo, show_progress = false)
        @unpack bas = data
        nlost = count(==(-1), bas)
        nlost > 0 && return NaN, NaN, NaN, file, "$nlost initial conditions labelled -1"
        Sb, Sbb = basin_entropy(bas, 10)
        return Sb, Sbb, length(unique(bas)), file, ""
    catch e
        e isa InterruptException && rethrow()
        return NaN, NaN, NaN, file, replace(sprint(showerror, e), r"\s+" => " ")
    end
end

# Compute basin_stats on grid g for every case, in parallel
function run_pass(cases, g, desc)
    Sb = fill(NaN, Np, Np)
    Sbb = fill(NaN, Np, Np)
    Na = fill(NaN, Np, Np)
    failed = Tuple{String, Float64, Float64, String}[]
    lk = ReentrantLock()
    prog = Progress(length(cases); desc)
    Threads.@threads :dynamic for c in cases
        j, k = Tuple(c)
        Sb[j,k], Sbb[j,k], Na[j,k], file, reason = basin_stats(j, k, g)
        isempty(reason) || lock(() -> push!(failed, (file, δrange[j], σrange[k], reason)), lk)
        next!(prog)
    end
    finish!(prog)
    return Sb, Sbb, Na, failed
end

# First pass: count the attractors on the coarse grid
all_cases = vec(CartesianIndices((Np, Np)))
_, _, Na_coarse, failed_coarse = run_pass(all_cases, grid_coarse, "Coarse pass: ")

# Second pass: fine basins where there are at least two attractors. Points that
# failed on the coarse grid (Na_coarse = NaN) are also recomputed.
fine_cases = filter(c -> !(Na_coarse[c] < 2), all_cases)
println("$(length(fine_cases)) of $(length(all_cases)) points have at least two attractors ",
    "or failed in the coarse pass ($(length(failed_coarse)) failed)")
Sb, Sbb, Na, failed = run_pass(fine_cases, grid, "Fine pass: ")

single = Na_coarse .== 1
Sb[single] .= 0.0
Sbb[single] .= 0.0
Na[single] .= 1.0

open(failed_file, "w") do io
    println(io, "file\tdelta\tsigma\treason")
    for (file, δ, σ, reason) in sort(failed)
        println(io, file, '\t', δ, '\t', σ, '\t', reason)
    end
end
println("$(length(failed)) of $(length(fine_cases)) fine basins failed, listed in $failed_file")

wsave(datadir("fig2_data.jld2"), @strdict Sb Sbb Na Na_coarse δrange σrange dps grid grid_coarse)
