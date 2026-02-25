close all
syms x
wAT1 = x;
wAT2 = x^2;
wWu = 2*x - x^2;
gAT = (1-x)^2;

C11 = 1;

gWu = ((1-x)^2)./(((1-x)^2)+3*x*(1-0.5*x));
gRat = (1-x)/((1-x)+1.5*x);

load('HexagonBenchmark03.mat')
gHexa = degradation.fun{1,1,1,1}(x);
load('HoneycombBenchmark03.mat')
gHoney = degradation.fun{1,1,1,1}(x);

figure(100)
% epsAT1 = sqrt(2*5e-3/((8/3)*0.1*C11))*sqrt(-diff(wAT1)/diff(gAT));
% fplot(epsAT1,[0 1])
% epsAT2 = sqrt(2*5e-3/(2*0.1*C11))*sqrt(-diff(wAT2)/diff(gAT));
% fplot(epsAT2,[0 1])
% epsWu =  sqrt(2*5e-3/(pi*0.1*C11))*sqrt(-diff(wWu)/diff(gWu));
% fplot(epsWu,[0 1])
% epsWuAT1 =  sqrt(2*5e-3/(pi*0.1*C11))*sqrt(-diff(wAT1)/diff(gWu));
% fplot(epsWuAT1,[0 1])
% epsRat = sqrt(2*5e-3/(pi*0.1*C11))*sqrt(-diff(wAT1)/diff(gRat));
% fplot(epsRat,[0 1])
epsHexa = sqrt(-diff(wAT1)/diff(gHexa));
a = fplot(epsHexa,[0 1]);
plot(a.YData,a.XData)

load('1ElemHexagonSimplified.mat');
hold on
plot(outputData.displacement.value,outputData.damage.maxValue)

% epsHoney = sqrt(2*5e-3/((8/3)*0.1*C11))*sqrt(-diff(wAT1)/diff(gHoney));
% fplot(epsHoney,[0 1])
legend('AT1','AT2','Wu','WuAT1','Rational','Hexagon','Honeycomb')

figure(101)
sigHexa =  sqrt(-gHexa*gHexa*diff(wAT1)/diff(gHexa));
fplot(sigHexa,[0 1]);

%% Computation 2D
syms x 
[dataHexa]  = load('HexagonBenchmark03.mat');
C11hexa = dataHexa.degradation.fun{1,1,1,1}(x);
C12hexa = dataHexa.degradation.fun{2,2,1,1}(x);
C33hexa = dataHexa.degradation.fun{1,2,1,2}(x);

bulkHexa   = (C11hexa-C33hexa);
shearHexa  = C33hexa;

[dataHoney] = load('HoneycombBenchmark03.mat');
C11honey = dataHoney.degradation.fun{1,1,1,1}(x);
C12honey = dataHoney.degradation.fun{2,2,1,1}(x);
C33honey = dataHoney.degradation.fun{1,2,1,2}(x);
bulkHoney  = (C11honey-C33honey);
shearHoney = C33honey;

k = bulkHexa;
m = shearHexa;
CPrime1Dhexa = ( (diff(k)+diff(m))*(1+((k-m)/(k+m))^2) - 2*(diff(k)-diff(m))*((k-m)/(k+m)));

k = bulkHoney;
m = shearHoney;
CPrime1Dhoney = ( (diff(k)+diff(m))*(1+((k-m)/(k+m))^2) - 2*(diff(k)-diff(m))*((k-m)/(k+m)));


epsHexa = sqrt(-diff(wAT1)/CPrime1Dhexa);
epsHoney = sqrt(-diff(wAT1)/CPrime1Dhoney);
figure(102)
a = fplot(epsHexa,[0 1]);
aXData = a.XData; aYData = a.YData;
b = fplot(epsHoney,[0 1]);
bXData = b.XData; bYData = b.YData;
plot(aYData,aXData)
hold on
plot(bYData,bXData)
legend('Hexa','Honey')
title('Damage-displacement')

figure(103)
hold on
fplot(CPrime1Dhexa,[0 1])
fplot(CPrime1Dhoney,[0 1])
legend('Hexa','Honey')
title('Degradation derivative equivalent to 1D')


