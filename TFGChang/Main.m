%% Test for TestNaca
tic; % Inicia el contador de tiempo
Naca.M     = 0.07;
Naca.p     = 0.8;
Naca.t     = 0.12;
Naca.chord = 2;
Naca.AoA   = 0;

NacaClass = TestNaca(Naca);
NacaClass.compute();
%NacaClass.validate();
%NacaClass.print();

%% Control Parameters.

m     = 0:0.005:0.09;
p     = 0.2:0.05:0.8;
t     = 0.05:0.05:0.4;
alpha = 0:1:30;

%% Code to compute Efficiency dataset

for k = 1:length(t)

    for i = 1:length(m)
    
        for j = 1:length(p)
        
            Naca.M     = m(i);
            Naca.p     = p(j);
            Naca.t     = t(k);
            Naca.chord = 1;
            Naca.AoA   = 10;
            
            NacaClass  = TestNaca(Naca);
            NacaClass.compute();
            %NacaClass.validate();
            NacaClass.print();
        
        end
    
    end
end

%% Plot efficiency vs p and m

data = load('EPerAoA0.txt');

E = data(:, end);
E = reshape(E, length(p), [])';

[P,M] = meshgrid(p,m);

surface(M,P,E);
shading interp;               
xlabel('Max Camber (m)');          
ylabel('Max Camber Position (p)');                 
title('L/D ratio versus Max Camber (m) and Max Camber Position (p)');
axis tight; 
colorbar;           
caxis([-0.08, 0.08]); 

%% Plot efficiency vs p and m

dataWPF = load('EPerAoA0WPF.txt');
data    = load('EPerAoA0.txt');

E = data(:, end);
E = reshape(E, length(p), [])';

EWPF = dataWPF(:, end);
EWPF = reshape(EWPF, length(p), [])';

ReError = abs(EWPF - E) / E;

[P,M] = meshgrid(p,m);

surface(M,P,ReError);
shading interp;               
xlabel('Max Camber (m)');          
ylabel('Max Camber Position (p)');                 
title('Relative Error of the L/D Ratio with and without a Pressure Filter as a Function of Maximum Camber (m) and Maximum Camber Position (p).');
axis tight; 
colorbar;           
caxis([-0.08, 0.08]); 


%% Code to compute L,D vs alpha dataset

for i = 1:length(alpha)

Naca.M     = 0.02;
Naca.p     = 0.4;
Naca.t     = 12/100;
Naca.chord = 1;
Naca.AoA   = alpha(i);

NacaClass  = TestNaca(Naca);
NacaClass.compute();
%NacaClass.validate();
NacaClass.print();
    
end

%% Plot L vs alpha and D vs L

data = load('results.txt');

L = data(:, end-2);
D = data(:, end-1);

figure;
plot(alpha, L);
xlabel('\alpha (deg)');
ylabel('L');
title('Lift vs AoA(\alpha)');
grid on;

figure;
plot(L, D);
xlabel('L');
ylabel('D');
title('Drag vs Lift');
grid on;   