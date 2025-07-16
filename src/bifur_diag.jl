function _get_basin_entropy(dps, δ, grid;  show_progress = false) 
        dps.δ = δ
        mapper = get_mapper(dps)
        bas, _ = basins_of_attraction(mapper, grid; show_progress) 
        Sb,Sbb = basin_entropy(bas, 10) 
        return Sb, Sbb, bas
end

function compute_entropy(d)
    @unpack  δrange, dps, grid = d
    Sb = zeros(length(δrange))
    Sbb = zeros(length(δrange))
    @Threads.threads for k in eachindex(δrange)
        @show δrange[k]
        Sb[k], Sbb[k], _ = _get_basin_entropy(deepcopy(dps), δrange[k], grid)
    end
    return @strdict(Sb, Sbb, δrange)
end

function get_entropy(δrange, dps, grid; force = false)
    δi = δrange[1]
    δf = δrange[end]
    N = length(δrange)
    d = @strdict  dps grid δrange δi δf N
    data, file = produce_or_load(
        compute_entropy,
        d, 
        datadir();
        prefix = "model_entropy", storepatch = false,
        suffix = "jld2", force = force
    )
    return data
end



function get_branches(prange, pidx, dps ,att)
    s = Vector{Int32}[]
    for e in att
        push!(s, collect(keys(e)))
    end
    s = unique(vcat(s...))
    branches = Dict( s[k] => Vector{Vector{Float64}}() for k in 1:length(s))

    # ptlst = Vector{Vector{Float64}}()
    smap = get_smap(dps)
    T = 2000
    Ttr = 500;
    for (k,el) in enumerate(att)
            for p in el
                set_parameters!(smap, [pidx => prange[k]])
                tra,t = trajectory(smap, T, p[2][1]; Ttr)
                  for y in tra[1900:2000]
                     v = [prange[k], y[1], y[2]]
                     push!(branches[p[1]], v)
                  end
            end
    end
    branches = Dict(k => StateSpaceSet(branches[k]) for k in keys(branches))
    return branches
end


function get_bif_points(branches::Dict)
    s = collect(keys(branches))
    bif_points = StateSpaceSet(branches[s[1]])
    for k in 2:length(s)
        append!(bif_points, StateSpaceSet(branches[s[k]]))
    end
    return bif_points
end

using Random: shuffle!, Xoshiro
function colors_from_keys(ukeys)
    # Unfortunately, until `to_color` works with `Cycled`,
    # we need to explicitly add here some default colors...
    COLORS = [
        "#7143E0",
        "#191E44",
        "#0A9A84",
        "#AF9327",
        "#5F166D",
        "#6C768C",
    ]
    if length(ukeys) ≤ length(COLORS)
        colors = [COLORS[i] for i in eachindex(ukeys)]
    else # keep colorscheme, but add extra random colors
        n = length(ukeys) - length(COLORS)
        colors = shuffle!(Xoshiro(123), collect(cgrad(COLORS, n+1; categorical = true)))
        colors = append!(to_color.(COLORS), colors[1:(end-1)])
    end
    return Dict(k => colors[i] for (i, k) in enumerate(ukeys))
end
