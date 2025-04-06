%% Compute Datas for Mesh Convergence Study

Naca.length = 8;
Naca.height = 4;

Naca.M      = 0.09;
Naca.p      = 0.8;
Naca.t      = 0.10;
Naca.chord  = 1;
Naca.AoA    = 10;

nx = 50:50:700;
E  = zeros(size(nx));

for i = 1:length(nx)

    disp(i);
    Naca.nx    = nx(i);

    NacaClass  = TestNaca(Naca);
    NacaClass.compute();
    
    E(i)       = NacaClass.E;     

end

save("data10.mat","E")
