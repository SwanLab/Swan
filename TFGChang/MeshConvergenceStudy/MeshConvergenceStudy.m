%% Compute Datas for Stokes Flow Mesh Convergence Study

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


%% Compute Datas for Navier-Stokes Flow Mesh Convergence Study

Naca.flowType  = "NavierStokes";
Naca.length    = 8;
Naca.height    = 4;
Naca.M         = 0.09;
Naca.p         = 0.8;
Naca.t         = 0.10;
Naca.chord     = 1;
Naca.AoA       = 10;
Naca.nu        = 0.1;


nx = 20:20:500;
E  = zeros(size(nx));

for i = 24:length(nx)

    disp(i);
    Naca.convectVel = 0;
    Naca.nx    = nx(i);

    NacaClass  = TestNaca(Naca);
    NacaClass.compute();
    NacaClass.print();

    E(i)       = NacaClass.E;

end

save("data10.mat","E")
