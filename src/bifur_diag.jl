

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
