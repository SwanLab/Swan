%% Compute Datas

Naca.length = 8;
Naca.height = 4;

Naca.M      = 0.04;
Naca.p      = 0.4;
Naca.t      = 0.12;
Naca.chord  = 1;
Naca.AoA    = 0;

nx = 50:50:800;
E  = zeros(size(nx));

for i = 1:length(nx)

    disp(i);
    Naca.nx    = nx(i);

    NacaClass  = TestNaca(Naca);
    NacaClass.compute();
    
    E(i)       = NacaClass.E;     

end

errorRelativeE = (E(:) - E(end)) ./ E(end);

%% Plot relative error vs nx

figure;
plot(nx, errorRelativeE);
grid on;
xlabel('Number of Elements in x-direction, n_x');
ylabel('Relative Efficiency Error');
title('Relative Efficiency Error vs. Number of Elements in x-direction');