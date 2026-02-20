syms x
wAT1 = x;
wAT2 = x^2;
wWu = 2*x - x^2;
gAT = (1-x)^2;

gWu = ((1-x)^2)./(((1-x)^2)+3*x*(1-0.5*x));
gRat = (1-x)/((1-x)+1.5*x);

figure(100)
epsAT1 = sqrt(5e-3/((8/3)*0.1))*sqrt(-diff(wAT1)/diff(gAT));
fplot(epsAT1,[0 1])

hold on

epsAT2 = sqrt(5e-3/(2*0.1))*sqrt(-diff(wAT2)/diff(gAT));
fplot(epsAT2,[0 1])

epsWu =  sqrt(5e-3/(pi*0.1))*sqrt(-diff(wWu)/diff(gWu));
fplot(epsWu,[0 1])

epsWuAT1 =  sqrt(5e-3/(pi*0.1))*sqrt(-diff(wAT1)/diff(gWu));
fplot(epsWuAT1,[0 1])

epsRat = sqrt(5e-3/(pi*0.1))*sqrt(-diff(wAT1)/diff(gRat));
fplot(epsRat,[0 1])

load('HexagonBenchmark03.mat')
gHexa = degradation.fun{1,1,1,1}(x);
epsHexa = sqrt(5e-3/(pi*0.1))*sqrt(-diff(wAT1)/diff(gHexa));
fplot(epsHexa,[0 1])

load('HoneycombBenchmark03.mat')
gHoney = degradation.fun{1,1,1,1}(x);
epsHoney = sqrt(5e-3/(pi*0.1))*sqrt(-diff(wAT1)/diff(gHoney));
fplot(epsHoney,[0 1])

legend('AT1','AT2','Wu','WuAT1','Rational','Hexagon','Honeycomb')