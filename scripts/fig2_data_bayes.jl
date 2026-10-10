# Data for Fig. 2 with the Bayesian update sampler: number of attractors, basin
# entropy and boundary basin entropy on a (δ, σ) grid. For each δ, attractors
# are continued along σ with AttractorSeedContinueMatch; the continuations for
# different δ run in parallel. Run with: julia -t auto fig2_data_bayes.jl
#
# The region [-5,5]² is tiled into n_tiles² boxes. At the first σ every box is
# sampled with dense_n initial conditions; afterwards with sparse_n, and a box
# is re-sampled densely only when its Bayes factor says the new labels do not
# fit its history (see BayesianUpdateSampler). Sb is the mean over boxes of the
# posterior mean entropy of each box, with its posterior variance var_Sb. Sbb
# is the same mean restricted to boundary boxes. Na counts every attractor found
# at each σ, including those found only by seeding from the previous σ.
#
# Each δ continuation is cached, so an interrupted run resumes. A continuation
# that throws, or a σ where some samples are labelled -1 (the mapper did not
# converge), gives NaN and is listed in data/fig2_bayes_failed.tsv. To redo them:
# tighten the tolerances in get_mapper / get_smap, set redo_failed = true and run
# again. Only the listed δ continuations are recomputed.
using DrWatson
@quickactivate
using Attractors
using ProgressMeter
using SpecialFunctions
using Statistics

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))
include(srcdir("inference_stuff.jl"))

# Boundary basin entropy: mean posterior entropy over the boxes where more than
# one basin has weight. A label keeps a decaying weight λⁿα after it was last
# seen in a box, so only labels with α ≥ αmin (about one recent observation) count.
function boundary_entropy(alphas; αmin = 1.0)
    S = [bayes_entropy(a) for a in alphas if count(≥(αmin), values(a)) > 1]
    return isempty(S) ? 0.0 : mean(S)
end

# Workaround for StateSpaceSets: statespace_sampler allocates Threads.nthreads()
# scratch buffers but indexes them with Threads.threadid(), which goes past
# nthreads() when Julia also runs an interactive thread (default with -t N on 1.12+)
function pad_generator_buffers!(sampler)
    for g in sampler.generators, _ in (length(g.dummies) + 1):Threads.maxthreadid()
        push!(g.dummies, similar(first(g.dummies)))
    end
    return sampler
end

function compute_sigma_continuation_bayes(d)
    @unpack dps, σrange, region, n_tiles, sparse_n, dense_n, λ, βprior, global_reset, seed = d

    # compute global continuation of attractors
    sampler = BayesianUpdateSampler(region, n_tiles;
        sparse_n, dense_n, λ, β = βprior, global_reset, seed, history = true)
    pad_generator_buffers!(sampler)
    mapper = get_mapper(dps)
    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pcurve = [Dict(:σ => σ) for σ in σrange]
    out = global_continuation(ascm, pcurve, sampler; show_progress = false)

    hist = bayesian_sampler_history(sampler)
    @assert length(hist.alphas) == length(σrange)
    est = bayes_estimates(hist)
    fractions_cont = out.fractions
    attractors_cont = out.attractors

    N = length(σrange)
    Sb = fill(NaN, N); var_Sb = fill(NaN, N); Sbb = fill(NaN, N); Na = fill(NaN, N)
    lost = falses(N)
    for i in 1:N
        lost[i] = haskey(fractions_cont[i], -1)
        lost[i] && continue
        Na[i] = length(attractors_cont[i])
        Sb[i] = est.mean_S[i]
        var_Sb[i] = est.var_S[i]
        Sbb[i] = boundary_entropy(hist.alphas[i])
    end
    n_panics = est.n_panics
    global_resets = est.global_resets
    resamplings = out.other["resamplings"]
    return @strdict(Sb, var_Sb, Sbb, Na, lost, n_panics, global_resets, resamplings,
        fractions_cont, attractors_cont, σrange, dps)
end

redo_failed = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
Np = 200
δrange = range(-30, 30, length = Np)
# Continue from low to high dissipation: attractors with small basins found at
# low σ are followed through seeding until they disappear
σrange = range(0.1, 0.5, length = Np)

# Bayesian sampler. βprior is the Dirichlet pseudo-count (β is taken by the model)
region = ((-5.0, 5.0), (-5.0, 5.0))
n_tiles = 1        # 25² boxes of side 0.4
sparse_n = 100       # ics per box at each σ: 6250 in total
dense_n = 5000  # ics per box at the first σ and when a box raises an alarm
λ = 0.7             # forgetting factor of the priors
βprior = 0.5
global_reset = true   # re-learn every box when the set of attractors changes
seed = 1
failed_file = datadir("fig2_bayes_failed.tsv")

bayes_key = (region, n_tiles, sparse_n, dense_n, λ, βprior, global_reset, seed)
bayes_file(dp) = datadir("fig2_bayes_" * hash_name(dp, σrange, bayes_key...) * ".jld2")

to_redo = Set{String}()
if redo_failed && isfile(failed_file)
    to_redo = Set(first(split(l, '\t')) for l in readlines(failed_file)[2:end])
    println("Recomputing $(length(to_redo)) failed δ continuations")
end

Sb = fill(NaN, Np, Np)
var_Sb = fill(NaN, Np, Np)
Sbb = fill(NaN, Np, Np)
Na = fill(NaN, Np, Np)
n_panics = zeros(Int, Np, Np)
failed = Tuple{String, Float64, Float64, String}[]
lk = ReentrantLock()
prog = Progress(Np)
Threads.@threads :dynamic for j in eachindex(δrange)
    dp = deepcopy(dps)
    dp.δ = δrange[j]
    file = bayes_file(dp)
    d = Dict("dps" => dp, "σrange" => σrange, "region" => region, "n_tiles" => n_tiles,
        "sparse_n" => sparse_n, "dense_n" => dense_n, "λ" => λ, "βprior" => βprior,
        "global_reset" => global_reset, "seed" => seed)
    try
        data, _ = produce_or_load(compute_sigma_continuation_bayes, d, datadir();
            filename = _ -> hash_name(dp, σrange, bayes_key...),
            prefix = "fig2_bayes", storepatch = false, suffix = "jld2",
            force = file in to_redo, verbose = false)
        Sb[j,:] .= data["Sb"]
        var_Sb[j,:] .= data["var_Sb"]
        Sbb[j,:] .= data["Sbb"]
        Na[j,:] .= data["Na"]
        n_panics[j,:] .= data["n_panics"]
        for (k, l) in enumerate(data["lost"])
            l && lock(() -> push!(failed, (file, dp.δ, σrange[k], "samples labelled -1")), lk)
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

wsave(datadir("fig2_data_bayes.jld2"), @strdict(Sb, var_Sb, Sbb, Na, n_panics, δrange, σrange,
    dps, region, n_tiles, sparse_n, dense_n, λ, βprior, global_reset, seed))
