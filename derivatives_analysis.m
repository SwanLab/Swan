a1 = (4/pi)*(5e-3*210)/(1.984^2 *0.1);
a2 = (4/pi)*(5e-3*210)/(1^2 *0.1);
C11 = 210/((1+0.3)*(1-0.3));
syms x
ratFun1 = C11*(1-x)/((1-x)+a1*x/2);
ratFun2 = C11*(1-x)/((1-x)+a2*x/2);
ratWu = C11*(1-x)^2 / ((1-x)^2 + a1*x*(1-0.5*x));

load('HexagonBenchmark03.mat');
C11hexa = degradation.fun{2,2,2,2};
dC11hexa = degradation.dfun{2,2,2,2};
d2C11hexa = degradation.ddfun{2,2,2,2};

load('HoneycombBenchmark03.mat');
C11honey = degradation.fun{1,1,1,1};
dC11honey = degradation.dfun{1,1,1,1};
d2C11honey = degradation.ddfun{1,1,1,1};

tiledlayout(1,3)
nexttile
hold on
fplot(ratFun1,[0 1])
fplot(ratFun2,[0 1])
fplot(ratWu,[0 1])
fplot(@(x) 210*C11hexa(x),[0 1])
fplot(@(x) 210*C11honey(x),[0 1])
title('Function')

nexttile
hold on
fplot(diff(ratFun1),[0 1])
fplot(diff(ratFun2),[0 1])
fplot(diff(ratWu),[0 1])
fplot(@(x) 210*dC11hexa(x),[0 1])
fplot(@(x) 210*dC11honey(x),[0 1])
title('First derivative')

nexttile
hold on
fplot(diff(diff(ratFun1)),[0 1])
fplot(diff(diff(ratFun2)),[0 1])
fplot(diff(diff(ratWu)),[0 1])
fplot(@(x) 210*d2C11hexa(x),[0 1])
fplot(@(x) 210*d2C11honey(x),[0 1])
title('Second derivative')

legend('Rational 1.98','Rational 1','Rational Wu 1.98','Hexagon','Honeycomb')
