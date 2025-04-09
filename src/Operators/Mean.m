function fM = Mean(f,quad)
Vol  = f.mesh.computeVolume();
intF = Integrator.compute(f,f.mesh,quad);
fM   = intF/Vol;
end