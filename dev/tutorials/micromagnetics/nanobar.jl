using JuMag
using Random
using NPZ

mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);

sim = Sim(mesh, driver="SD", name="bar");

set_Ms(sim, 8e5)

add_exch(sim, 1.3e-12, name="exch")
add_demag(sim, name="demag");

Random.seed!(12345)
init_m0_random(sim);

relax(sim, maxsteps=2000, stopping_dmdt=0.01)

npzwrite("bar.npy", Array(sim.spin))

save_vtk(sim, "bar", fields=["exch", "demag"])

using CairoMakie

function plot_spatial_m()
    folder = @__DIR__
    m = npzread(folder*"/bar.npy")

    nx, ny, nz = 30, 5, 2
    xs = [i*2 for i=1:nx]
    ys = [j*2 for j=1:ny]
    m = reshape(m, 3, nx, ny, nz)
    mx = m[1,:,:,1]
    my = m[2,:,:,1]

    fig = Figure(resolution = (800, 200))
    ax = Axis(fig[1, 1], backgroundcolor = "white")

    arrows!(ax, xs, ys, mx, my, arrowsize = 10, lengthscale = 2,
            arrowcolor = vec(my), linecolor = vec(my), align = :center)


    return fig

end

plot_spatial_m()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

