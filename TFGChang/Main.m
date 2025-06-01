%% Test for TestNaca
Naca.flowType  = "Stokes";
Naca.length      = 8;
Naca.height      = 4;
Naca.nx          = 420;
Naca.M           = 0.02;
Naca.p           = 0.4;
Naca.t           = 0.12;
Naca.chord       = 1;
Naca.AoA         = 10;

NacaClass = TestNaca(Naca);
NacaClass.compute();
%NacaClass.validate();
%NacaClass.print();

%system('shutdown /s /t 60');

%% Test for TestNaca NV
tic
nuRef = [1, 0.1, 1/20, 1/30,1/40, 1/50];
Naca.convectVel = 0;

for i = 1 : length(nuRef)
    %length(nuRef)
tic
Naca.flowType  = "NavierStokes";
Naca.length    = 8;
Naca.height    = 4;
Naca.nx        = 210;
Naca.M         = 0.02;
Naca.p         = 0.4;
Naca.t         = 0.12;
Naca.chord     = 1;
Naca.AoA       = 5;
Naca.nu        = nuRef(i);
%Naca.i         = i;
%0.755;
disp(i)
disp(Naca.nu)

NacaClass = TestNaca(Naca);
NacaClass.compute();
%NacaClass.validate();
NacaClass.print();

Naca.convectVel = NacaClass.velocityFun.fValues;

end
toc
%% Control Parameters.

Naca.length = 8;
Naca.height = 4;
Naca.nx     = 420;
m     = 0.01:0.01:0.09;
p     = 0.2:0.1:0.8;
t     = 0.1:0.02:0.4;
alpha = -15:1:15;

%% Code to compute Efficiency dataset 

% for Symetrical Airfoil

for k = 1:length(t)

            fprintf("m = %f\n", 0.0);
            fprintf("p = %f\n", 0.2);
            fprintf("t = %f\n", t(k));
            Naca.M     = 0.0;
            Naca.p     = 0.2;
            Naca.t     = t(k);
            %t(k);
            Naca.chord = 1;
            Naca.AoA   = 12.5;

            NacaClass  = TestNaca(Naca);
            NacaClass.compute();
            NacaClass.print();
end

% for Asymetrical Airfoil

for k = 1:length(t)

    for i = 1:length(m)
    
        for j = 1:length(p)
        
            fprintf("m = %f\n", m(i));
            fprintf("p = %f\n", p(j));
            fprintf("t = %f\n", t(k));
            Naca.M     = m(i);
            Naca.p     = p(j);
            Naca.t     = t(k);
            Naca.chord = 1;
            Naca.AoA   = 12.5;
            
            NacaClass  = TestNaca(Naca);
            NacaClass.compute();
            NacaClass.print();
            
        end
    
    end
end


%% Code to compute L,D vs alpha dataset
alpha = 10:0.5:12.5;

for i = 1:length(alpha)

fprintf("alpha = %f\n", alpha(i));
Naca.flowType  = "Stokes";
Naca.length    = 8;
Naca.height    = 4;
Naca.nx        = 420;
Naca.M         = 0;
Naca.p         = 0;
Naca.t         = 0.1;
Naca.chord     = 1;
Naca.AoA   = alpha(i);

NacaClass  = TestNaca(Naca);
NacaClass.compute();
%NacaClass.validate();
NacaClass.print();
    
end

%% Plot efficiency vs p and m

data = load('E_AOA0_nx600_PF_FL.txt');

E = data(:, end);
E = reshape(E, length(p), [])';

[P,M] = meshgrid(p,m);

figure;
surface(M,P,E);
shading interp;               
xlabel('Max Camber (m)');          
ylabel('Max Camber Position (p)');                 
title('L/D ratio versus Max Camber (m) and Max Camber Position (p)');
axis tight; 
colorbar;           
caxis([-0.08, 0.08]); 
hold on;
mesh(M, P, E, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.5); % Cuadrícula negra

%% Plot Relative Error of the L/D Ratio with and without a Pressure Filter

dataWPF = load('E_AOA0_nx600_PF_FL.txt');
data    = load('E_AOA0_nx600_PNF_FL.txt');

E = data(:, end);
E = reshape(E, length(p), [])';

EWPF = dataWPF(:, end);
EWPF = reshape(EWPF, length(p), [])';

ReError = abs(EWPF - E) ./ E;

[P,M] = meshgrid(p,m);

figure;
surface(M,P,ReError);
shading interp;               
xlabel('Max Camber (m)');          
ylabel('Max Camber Position (p)');                 
title('Relative Error of the L/D Ratio with and without a Pressure Filter as a Function of Maximum Camber (m) and Maximum Camber Position (p).');
axis tight; 
colorbar;           
%caxis([-0.08, 0.08]); 

%% Plot L vs alpha and D vs L

data = load('LD_N2412_AoA-+15_nx420_PF_NFL.txt');

L = data(:, end-2);
D = data(:, end-1);
E = data(:,end);

% Lift figure
figure;
plot(alpha, L, 'b');
xlabel('\alpha (deg)');
ylabel('Lift (L)');
title('Lift vs AoA (\alpha)');
xlim([-15 15]);
grid on;

% Drag figure
figure;
plot(alpha, D, 'r');
xlabel('\alpha (deg)');
ylabel('Drag (D)');
title('Drag vs AoA (\alpha)');
xlim([-15 15]);
grid on;

% Efficiency figure
figure;
plot(alpha, E, 'g');
xlabel('\alpha (deg)');
ylabel('Efficiency (E)');
title('Efficiency vs AoA (\alpha)');
xlim([-15 15]);
grid on;

%% Plot L vs alpha and D vs L (Pressure Filtered and No Filtered)

dataPF = load('LD_N2412_AOA-+15_nx600_PF_FL.txt');
dataPNF = load('LD_N2412_AOA-+15_nx600_PNF_FL.txt');

LPF = dataPF(:, end-2);
DPF = dataPF(:, end-1);

LPNF = dataPNF(:, end-2);
DPNF = dataPNF(:, end-1);

figure;
plot(alpha, LPF);
hold on;
plot(alpha, LPNF);
legend('Filtered',"No Filtered", 'Location', 'best')
xlabel('\alpha (deg)');
ylabel('L');
title('Lift vs AoA(\alpha)');
grid on;

figure;
plot(LPF, DPF);
hold on;
plot(LPNF, DPNF);
legend('Filtered',"No Filtered", 'Location', 'best')
xlabel('L');
ylabel('D');
title('Drag vs Lift');
grid on;   

