using JuMag
using Printf
using NPZ

JuMag.cuda_using_double(true)

function m0_fun(i, j, k, dx, dy, dz)
    r = 25
    if ((i - 80)^2 + (j - 40)^2 < r^2)
        return (0.05, 0.01, -1)
    end
    return (0,0,1)
end

function relax_system()
    J = 1*meV
    D = 0.18*J

    mesh = CubicMeshGPU(nx=166, ny=80, nz=1, dx=0.5e-9, dy=0.5e-9, pbc="xy")

    sim = Sim(mesh, driver="SD", name="skx")
    set_mu_s(sim, mu_s_1)

    add_exch(sim, J, name="exch")
    add_dmi(sim, D, name="dmi")

    Hz= 1.5e-2*J / mu_s_1
    add_zeeman(sim, (0,0,Hz))
    init_m0(sim, m0_fun)

    relax(sim, maxsteps=1000, stopping_dmdt=1e-5, using_time_factor=false)

    save_vtk(sim, "skx")

    Q = compute_skyrmion_number(Array(sim.spin),mesh)

    return Q
end

relax_system()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

