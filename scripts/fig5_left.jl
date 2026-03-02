using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using ChaosTools

include(srcdir("model_mapper.jl"))
include(srcdir("bifur_diag.jl"))
function compute_delta_sweep(params::Dict)
    @unpack  δrange, Np, Nsamples, dps, grid= params

    # compute globl continuation of attractors
    sampler, _ = statespace_sampler(grid) 
    mapper = get_mapper(dps)

    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pidx = :δ

    fractions_cont, attractors_cont = global_continuation(
        ascm, δrange, pidx, sampler; samples_per_parameter = Nsamples
    )

    branches = get_branches(δrange, pidx, dps, attractors_cont)

    return @strdict(fractions_cont,  attractors_cont, δrange, dps, branches)
end

force = false
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = 1.0; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
res = 400
yg = range(-5,5, length=res); grid = (yg, yg) 
Np = 800; Nsamples = 3000
δrange = range(-30, -20, length = Np)


# Compute the fractions and attractor branches
params = @strdict Np Nsamples δrange dps grid
dat, _ = produce_or_load(compute_delta_sweep, params, datadir(); prefix = "delta_sweep_zoom", force)

# Get variables from dictionary
@unpack fractions_cont, attractors_cont, branches = dat

labs_args = (ylabel = L"i", xlabel = L"q", yticklabelsize = 20,  ylabelsize = 25, xticklabelsize = 20,  xlabelsize = 25)

ind = 400
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = δrange[ind]; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
# Compute the fractions and attractor branches
fig = Figure(size = (1000,400))

ax = Axis(fig[1,1]; labs_args...  )
u0 = attractors_cont[ind][1][1] 
T = 100
y,t = get_trajectory_map(u0, T, dps)
y2,t = get_trajectory_ode(u0, T, dps)
lines!(ax, y2[:,1], y2[:,2]; linewidth = 0.5, color = :blue)
scatter!(ax, y[:,1],y[:,2]; markersize = 10.0, color = :red)

labs_iq = (xticklabelsvisible = true, yticklabelsvisible = true, xlabel = L"\tau", xlabelsize = 25, ylabel = L"q", ylabelsize = 25, yticklabelsize = 20, xticklabelsize = 20)

ax_inset = Axis(fig[2, 1];   xticks = ([0, 1], ["0", "1"]), labs_iq...)
lines!(ax_inset, t[1:1000]./(2*π), y2[1:1000,1]; linewidth = 0.5, color = :blue)

ax = Axis(fig[1,2]; labs_args...)
u0 = attractors_cont[ind][6][1] 
T = 100
y,t = get_trajectory_map(u0, T, dps)
y2,t = get_trajectory_ode(u0, T, dps)
lines!(ax, y2[:,1], y2[:,2]; linewidth = 0.5, color = :red)
scatter!(ax, y[:,1],y[:,2]; markersize = 10.0, color = :red)

ax_inset = Axis(fig[2, 2];   xticks = ([0, 1, 2, 3], ["0", "1", "2", "3"]), labs_iq...)
lines!(ax_inset, t[1:2500]./(2*π), y2[1:2500,1]; linewidth = 0.5, color = :red)

ax = Axis(fig[1,3]; labs_args...)
u0 = attractors_cont[ind][5][1] 
T = 100
y,t = get_trajectory_map(u0, T, dps)
y2,t = get_trajectory_ode(u0, T, dps)
lines!(ax, y2[:,1], y2[:,2]; linewidth = 0.5, color = :black)
scatter!(ax, y[:,1],y[:,2]; markersize = 10.0, color = :red)

ax_inset = Axis(fig[2, 3]; labs_iq...)
lines!(ax_inset, t[1:10000]./(2*π), y2[1:10000,1]; linewidth = 0.5, color = :black)

# save("orbit_delta=-25.00.png",fig)
save("fig5a_1.png",fig)


ind = 325
σ = 0.3; ω = 1.0; μ = 35.0; η = 0.08; δ = δrange[ind]; β = 0.4
dps = model_parameters(ω, σ, β, η, μ, δ)
# Compute the fractions and attractor branches
fig = Figure(size = (1000,400))

ax = Axis(fig[1,1]; labs_args...)
u0 = attractors_cont[ind][1][1] 
T = 100
y,t = get_trajectory_map(u0, T, dps)
y2,t = get_trajectory_ode(u0, T, dps)
lines!(ax, y2[:,1], y2[:,2]; linewidth = 0.5, color = :blue)
scatter!(ax, y[:,1],y[:,2]; markersize = 10.0, color = :red)

ax_inset = Axis(fig[2, 1];   xticks = ([0, 1], ["0", "1"]), labs_iq...)
lines!(ax_inset, t[1:1000]./(2*π), y2[1:1000,1]; linewidth = 0.5, color = :blue)

ax = Axis(fig[1,2]; labs_args...)
u0 = attractors_cont[ind][4][1] 
T = 100
y,t = get_trajectory_map(u0, T, dps)
y2,t = get_trajectory_ode(u0, T, dps)
lines!(ax, y2[:,1], y2[:,2]; linewidth = 0.5, color = :red)
scatter!(ax, y[:,1],y[:,2]; markersize = 10.0, color = :red)

ax_inset = Axis(fig[2, 2];   xticks = ([0, 1], ["0", "1"]), labs_iq...)
lines!(ax_inset, t[1:1000]./(2*π), y2[1:1000,1]; linewidth = 0.5, color = :red)

ax = Axis(fig[1,3]; labs_args...)
u0 = attractors_cont[ind][5][1] 
T = 100
y,t = get_trajectory_map(u0, T, dps)
y2,t = get_trajectory_ode(u0, T, dps)
lines!(ax, y2[:,1], y2[:,2]; linewidth = 0.5, color = :black)
scatter!(ax, y[:,1],y[:,2]; markersize = 10.0, color = :red)

ax_inset = Axis(fig[2, 3];  labs_iq...)
lines!(ax_inset, t[1:1500]./(2*π), y2[1:1500,1]; linewidth = 0.5, color = :black)

# save("orbit_delta=-25.94.png",fig)
save("fig5b_2.png",fig)

