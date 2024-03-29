function save_ovf(sim::AbstractSimGPU, fname::String; 
  type::DataType = _cuda_using_double.x ? Float64 : Float32) 

  T = _cuda_using_double.x ? Float64 : Float32
  mesh = sim.mesh
  nx,ny,nz = mesh.nx,mesh.ny,mesh.nz
  nxyz = nx*ny*nz

  ovf = OVF2{T}()
  ovf.xnodes = mesh.nx
  ovf.ynodes = mesh.ny
  ovf.znodes = mesh.nz
  ovf.xstepsize = mesh.dx
  ovf.ystepsize = mesh.dy
  ovf.zstepsize = mesh.dz
  ovf.type = type
  ovf.name = fname
  ovf.data = zeros(T,3*nxyz)
  copyto!(ovf.data, sim.spin)

  save_ovf(ovf, fname)
end


function read_ovf(sim::AbstractSimGPU, fname::String)
  T = _cuda_using_double.x ? Float64 : Float32
  ovf  = read_ovf(fname, T=T)
  nxyz = ovf.xnodes*ovf.ynodes*ovf.znodes
  if nxyz != sim.nxyz
      error("The ovf does not match sim.mesh")
  end
  copyto!(sim.prespin, ovf.data)
  copyto!(sim.spin, ovf.data)
end