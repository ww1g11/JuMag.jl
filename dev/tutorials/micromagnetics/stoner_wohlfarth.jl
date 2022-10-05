using JuMag

function relax_sw_loop()
    #We create a mesh for a cubic geometry 4nm x 4nm x 4nm
    mesh = FDMeshGPU(nx=4, ny=4, nz=4, dx=1e-9, dy=1e-9, dz=1e-9)

    #We create a simulation with the SD driver.
    sim = Sim(mesh, name="sw_cubic", driver="SD")

    #We set the saturation magnetization of the system.
    set_Ms(sim, 1.0e6)

    #Set the exchange constant
    add_exch(sim, 1.3e-11)

    #Set the initial state of the system.
    init_m0(sim, (-1, 1, 0))

    #Add the anisotropy
    add_anis(sim, 5e4, axis=(1,1,0))

    #Set the zeeman field.
    add_zeeman(sim,(0,0,0))

    #For each field, we relax the system to obtain its equilibrium state.
    for i=-100:5:100
        Hx = i*mT # A/m

        update_zeeman(sim, (Hx,0,0))

        #Relax the system with stopping_dmdt=0.05, the write_data function will be called if save_m_every is positive
        relax(sim, maxsteps=10000, stopping_dmdt=0.05, save_m_every=10000)

    end
end

relax_sw_loop()

using DelimitedFiles
using CairoMakie

function plot_loop()
    folder = @__DIR__
    data = readdlm(folder*"/sw_cubic_sd.txt", skipstart=2)
    m, H = data[:, 3], data[:, 8]

    fig = Figure(resolution = (800, 500))
    ax = Axis(fig[1, 1],
        xlabel = "H (A/m)",
        ylabel = "mx"
    )

    scatterlines!(ax, H, m, markersize = 6, color = :blue, markercolor = :orange)
    scatterlines!(ax, -H, -m, markersize = 6, color = :blue, markercolor = :orange)

    expected = 39788.736 # A/m
    vlines!(ax, [expected, -expected], color = :red, linestyle = :dash)

    return fig

end

plot_loop()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

