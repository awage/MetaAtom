# Data for Fig. 2 from sampled basins: number of attractors, basin entropy and
# boundary basin entropy on a (δ, σ) grid. For each δ, attractors are continued
# along σ with AttractorSeedContinueMatch; the continuations for different δ
# run in parallel. Run with: julia -t auto fig2_data_cont.jl
#
# The same Nsamples random initial conditions in [-5,5]² are used at every
# (δ, σ). Sb and Sbb are estimated from their labels with the nearest-neighbour
# estimator basin_entropy(points, labels, nneigh). Na counts every attractor
# found, including those found only by seeding from the previous σ, whose
# basins may be too small to be hit by the samples.
#
# Each δ continuation is cached, so an interrupted run resumes. A continuation
# that throws, or a σ where some samples are labelled -1 (the mapper did not
# converge), gives NaN and is listed in data/fig2_cont_failed.tsv. To redo them:
# tighten the tolerances in get_mapper / get_smap, set redo_failed = true and run
# again. Only the listed δ continuations are recomputed.
using DrWatson
@quickactivate
using Attractors
using ProgressMeter
using Random: Xoshiro

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))

# Progress log shared by all threads: a "start" line when the mapping at a
# (δ, σ) begins and a "done" line with its wall time when it ends. A "start"
# without "done" is a point still running; watch it with tail -f.
const progress_log = datadir("fig2_cont_progress.log")
const log_lock = ReentrantLock()
function log_progress(msg)
    lock(log_lock) do
        open(io -> println(io, Libc.strftime("%F %T", time()), "  t", Threads.threadid(), "  ", msg),
            progress_log, "a")
    end
end

# Sampler that always returns the same initial conditions and keeps their
# labels at each parameter value of the continuation. It also logs the time
# spent mapping the initial conditions at each σ.
struct LabelRecorder{V} <: InitialConditionsSampler
    ics::V
    labels::Vector{Vector{Int}}
    δ::Float64
    σrange::Vector{Float64}
    tstart::Base.RefValue{Float64}
end
LabelRecorder(ics, δ, σrange) = LabelRecorder(ics, Vector{Int}[], δ, collect(σrange), Ref(0.0))
current_σ(s::LabelRecorder) = s.σrange[length(s.labels) + 1]
function Attractors.generate_ics(s::LabelRecorder, args...)
    s.tstart[] = time()
    log_progress("start δ = $(s.δ), σ = $(current_σ(s))")
    return s.ics
end
Base.length(s::LabelRecorder) = length(s.ics)
function Attractors.update_sampler!(s::LabelRecorder, labels, args...)
    dt = time() - s.tstart[]
    log_progress("done  δ = $(s.δ), σ = $(current_σ(s)), $(round(dt; digits = 1)) s, " *
        "$(length(unique(labels))) labels, $(count(==(-1), labels)) lost")
    push!(s.labels, copy(labels))
end

function compute_sigma_continuation(d)
    @unpack dps, σrange, ics, nneigh = d

    # compute global continuation of attractors
    sampler = LabelRecorder(ics, dps.δ, σrange)
    mapper = get_mapper(dps)
    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pcurve = [Dict(:σ => σ) for σ in σrange]
    out = global_continuation(ascm, pcurve, sampler; show_progress = false)
    @assert length(sampler.labels) == length(σrange)

    N = length(σrange)
    Sb = fill(NaN, N); Sbb = fill(NaN, N); Na = fill(NaN, N)
    nlost = zeros(Int, N)
    for (i, labels) in enumerate(sampler.labels)
        nlost[i] = count(==(-1), labels)
        nlost[i] > 0 && continue
        Na[i] = length(out.attractors[i])
        Sb[i], Sbb[i] = basin_entropy(ics, labels, nneigh)
    end
    fractions_cont = out.fractions
    attractors_cont = out.attractors
    return @strdict(Sb, Sbb, Na, nlost, fractions_cont, attractors_cont, σrange, dps)
end

redo_failed = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
Np = 200
δrange = range(-30, 30, length = Np)
# Continue from low to high dissipation: attractors with small basins found at
# low σ are followed through seeding until they disappear
σrange = range(0.1, 0.5, length = Np)
Nsamples = 20000   # initial conditions per (δ, σ)
nneigh = 100       # neighbours in the entropy estimate
seed = 1
rng = Xoshiro(seed)
ics = StateSpaceSet([10 .* rand(rng, 2) .- 5 for _ in 1:Nsamples])
failed_file = datadir("fig2_cont_failed.tsv")

# ics is fully determined by (Nsamples, seed), so these go into the hash instead
cont_file(dp) = datadir("fig2_cont_" * hash_name(dp, σrange, Nsamples, seed, nneigh) * ".jld2")

to_redo = Set{String}()
if redo_failed && isfile(failed_file)
    to_redo = Set(first(split(l, '\t')) for l in readlines(failed_file)[2:end])
    println("Recomputing $(length(to_redo)) failed δ continuations")
end

Sb = fill(NaN, Np, Np)
Sbb = fill(NaN, Np, Np)
Na = fill(NaN, Np, Np)
failed = Tuple{String, Float64, Float64, String}[]
lk = ReentrantLock()
prog = Progress(Np)
Threads.@threads :dynamic for j in eachindex(δrange)
    dp = deepcopy(dps)
    dp.δ = δrange[j]
    file = cont_file(dp)
    d = Dict("dps" => dp, "σrange" => σrange, "ics" => ics, "nneigh" => nneigh)
    try
        data, _ = produce_or_load(compute_sigma_continuation, d, datadir();
            filename = _ -> hash_name(dp, σrange, Nsamples, seed, nneigh),
            prefix = "fig2_cont", storepatch = false, suffix = "jld2",
            force = file in to_redo, verbose = false)
        Sb[j,:] .= data["Sb"]
        Sbb[j,:] .= data["Sbb"]
        Na[j,:] .= data["Na"]
        for (k, n) in enumerate(data["nlost"])
            n > 0 && lock(() -> push!(failed, (file, dp.δ, σrange[k], "$n samples labelled -1")), lk)
        end
    catch e
        e isa InterruptException && rethrow()
        reason = replace(sprint(showerror, e), r"\s+" => " ")
        lock(() -> push!(failed, (file, dp.δ, NaN, reason)), lk)
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
nδ = length(unique(first.(failed)))
println("$(length(failed)) failures in $nδ of $Np δ continuations, listed in $failed_file")

wsave(datadir("fig2_data_cont.jld2"), @strdict Sb Sbb Na δrange σrange dps Nsamples nneigh seed)
