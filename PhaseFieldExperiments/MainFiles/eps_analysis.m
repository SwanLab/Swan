close all
syms x
wAT1 = x;
wAT2 = x^2;
wWu = 2*x - x^2;
gAT = (1-x)^2;

C11 = 1;

gWu = ((1-x)^2)./(((1-x)^2)+3*x*(1-0.5*x));
gRat = (1-x)/((1-x)+1.5*x);

load('HorizontalCrack.mat')
gHomog = degradation.fun{2,2,2,2}(x);

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
epsHomog = sqrt(-diff(wAT1)/diff(gHomog));
a = fplot(epsHomog,[0 1]);
plot(a.YData,a.XData)
legend('AT1','AT2','Wu','WuAT1','Rational','Homoggon','Honeycomb')

figure(101)
sigHomog =  sqrt(-gHomog*gHomog*diff(wAT1)/diff(gHomog));
fplot(sigHomog,[0 1]);

%% Computation 2D (Isotropic materials)
syms x 
[dataHomog]  = load('HorizontalCrack.mat');
C11Homog = dataHomog.degradation.fun{1,1,1,1}(x);
C12Homog = dataHomog.degradation.fun{2,2,1,1}(x);
C33Homog = dataHomog.degradation.fun{1,2,1,2}(x);

bulkHomog   = (C11Homog-C33Homog);
shearHomog  = C33Homog;


k = bulkHomog;
m = shearHomog;
CPrime1DHomog = ( (diff(k)+diff(m))*(1+((k-m)/(k+m))^2) - 2*(diff(k)-diff(m))*((k-m)/(k+m)));


epsHomog = sqrt(-diff(wAT1)/CPrime1DHomog);
figure(102)
a = fplot(epsHomog,[0 1]);
aXData = a.XData; aYData = a.YData;
plot(aYData,aXData)
legend('Homog')
title('Damage-displacement')

figure(103)
hold on
fplot(CPrime1DHomog,[0 1])
legend('Homog')
title('Degradation derivative equivalent to 1D')

%% Computation 2D (Anisotropic materials)


