using JuMag

mesh =  FDMesh(dx=5e-9, dy=5e-9, dz=5e-9, nx=1, ny=1, nz=1)

sim = Sim(mesh);

set_Ms(sim, 8.6e5)

init_m0(sim, (1,1,1))

demag = add_demag(sim);

zeeman = add_zeeman(sim, (0,0,1e5))

JuMag.effective_field(sim, sim.spin, 0.0)

println("Demag Energy: ",sum(demag.energy), " J")
println("Zeeman Energy: ",sum(zeeman.energy), " J")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
