using .CairoMakie
CairoMakie.activate!(type = "png")

function plot_poincare(bfield::MagneticField,
                       soln_ensemble;
                       resolution = (800, 800),
                       markersize = 2,
                      )
    trajectories = soln_ensemble.u
    T = eltype(soln_ensemble)
    r_vals = Vector{Vector{T}}(undef, length(trajectories))
    z_vals = Vector{Vector{T}}(undef, length(trajectories))

    for i in 1:length(trajectories)
        r_vals[i] = [ui[1] for ui in trajectories[i].u]
        z_vals[i] = [ui[2] for ui in trajectories[i].u]
    end

    r_max = maximum(maximum.(r_vals))
    r_min = minimum(minimum.(r_vals))
    r_buffer = 0.05*(r_max - r_min)
    z_max = maximum(maximum.(z_vals))
    z_min = minimum(minimum.(z_vals))
    z_buffer = 0.05*(r_max - r_min)
    fig = Figure(; resolution = resolution)
    ax = Axis(fig[1, 1], xlabel = "R (m)", ylabel = "Z (m)",
              limits = ((r_min - r_buffer, r_max + r_buffer),
              (z_min - z_buffer, z_max + z_buffer)),
              title = "Poincare")
    for i in eachindex(r_vals, z_vals)
        scatter!(ax, r_vals[i], z_vals[i], markersize = markersize)
    end
    fig
end
