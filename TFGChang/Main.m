%% Test for TestNaca

Naca.M     = 0.02;
Naca.p     = 0.4;
Naca.t     = 12/100;
Naca.chord = 1;
Naca.AoA   = 0;

NacaClass = TestNaca(Naca);
NacaClass.compute();
NacaClass.validate();
NacaClass.print();

%% Code to compute Efficiency dataset
m = 0:0.01:0.09;
p = 0.2:0.1:0.8;

for i = 1:length(m)

    for j = 1:length(p)
    
        Naca.M     = m(i);
        Naca.p     = p(j);
        Naca.t     = 12/100;
        Naca.chord = 1;
        Naca.AoA   = 0;
        
        NacaClass  = TestNaca(Naca);
        NacaClass.compute();
        %NacaClass.validate();
        NacaClass.print();
    
    end

end

%% Plot efficiency vs p and m

data = load('results.txt');

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

%% Code to compute L,D vs alpha dataset
alpha = 0:1:30;

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