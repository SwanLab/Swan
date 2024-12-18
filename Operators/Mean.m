function fM = Mean(f,mesh,quad)
Vol  = mesh.computeVolume();
intF = Integrator.compute(f,mesh,quad);
fM   = intF/Vol;
end