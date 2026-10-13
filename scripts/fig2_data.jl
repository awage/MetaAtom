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
failed_file = datadir("fig2_failed.tsv")

# to_redo = Set{String}()
# if redo_failed && isfile(failed_file)
#     to_redo = Set(first(split(l, '\t')) for l in readlines(failed_file)[2:end])
#     println("Recomputing $(length(to_redo)) failed basins")
# end

# Basin of the point (j, k) on grid g, and its Sb, Sbb, Na and attractors. A
# failure is returned as a non empty reason, with NaN values. att0 seeds the
# mapper with known attractors (see get_mapper).
function basin_stats(j, k, g; att0 = nothing)
    dp = deepcopy(dps)
    dp.δ = δrange[j]
    dp.σ = σrange[k]
    file = basins_file(dp, g)
    try
        data = get_basins(dp, g; force, show_progress = false, att0)
        @unpack bas, att = data
        nlost = count(==(-1), bas)
        nlost > 0 && return NaN, NaN, NaN, att, file, "$nlost initial conditions labelled -1"
        Sb, Sbb = basin_entropy(bas, 10)
        return Sb, Sbb, length(unique(bas)), att, file, ""
    catch e
        e isa InterruptException && rethrow()
        return NaN, NaN, NaN, nothing, file, replace(sprint(showerror, e), r"\s+" => " ")
    end
end

# Compute basin_stats on grid g for every case, in parallel. Seeds[j, k] are the
# attractors to seed the mapper with, or nothing. :greedy hands out one case at
# a time: the cost per case is very uneven (slow at low σ), and :dynamic would
# leave most threads idle at the end.
function run_pass(cases, g, desc; seeds = fill(nothing, Np, Np))
    Sb = fill(NaN, Np, Np)
    Sbb = fill(NaN, Np, Np)
    Na = fill(NaN, Np, Np)
    Att = Matrix{Any}(nothing, Np, Np)
    failed = Tuple{String, Float64, Float64, String}[]
    lk = ReentrantLock()
    prog = Progress(length(cases); desc)
    Threads.@threads :greedy for c in cases
        j, k = Tuple(c)
        Sb[j,k], Sbb[j,k], Na[j,k], Att[j,k], file, reason = basin_stats(j, k, g; att0 = seeds[j,k])
        isempty(reason) || lock(() -> push!(failed, (file, δrange[j], σrange[k], reason)), lk)
        next!(prog)
    end
    finish!(prog)
    return Sb, Sbb, Na, Att, failed
end

# First pass: count the attractors on the coarse grid
all_cases = vec(CartesianIndices((Np, Np)))
_, _, Na_coarse, Att_coarse, failed_coarse = run_pass(all_cases, grid_coarse, "Coarse pass: ")

# Second pass: fine basins where there are at least two attractors. Points that
# failed on the coarse grid (Na_coarse = NaN) are also recomputed.
fine_cases = filter(c -> !(Na_coarse[c] < 2), all_cases)
println("$(length(fine_cases)) of $(length(all_cases)) points have at least two attractors ",
    "or failed in the coarse pass ($(length(failed_coarse)) failed)")
# The fine pass is seeded with the attractors of the coarse pass (same mapper).
Sb, Sbb, Na, _, failed = run_pass(fine_cases, grid, "Fine pass: "; seeds = Att_coarse)

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

