using JuMag
using Test

include("test_mesh.jl")
include("test_zeeman.jl")
include("test_anis.jl")
include("test_exch.jl")
include("test_dmi.jl")
include("test_demag.jl")
include("test_fields.jl")
include("test_sim.jl")
include("test_llg.jl")
include("test_integrator.jl")
include("test_relax.jl")


include("atomistic/test_mesh.jl")
include("atomistic/test_fields.jl")
