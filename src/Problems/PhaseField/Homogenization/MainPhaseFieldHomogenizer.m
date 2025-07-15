clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0;
s.meshType   = 'Square';
s.meshN      = 100;
s.holeType   = 'Square';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

% [f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
% degradation.fun = f;
% degradation.dfun = df;
% degradation.ddfun = ddf;
% fplot(f,[0 1])
%save('CirclePerimeter','mat','phi','degradation')

%% DERIVATIVE
load('SquareAreaDerivativeNu0Sig1.mat')
[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('SquareAreaDerivativeNu0Sig1','mat','phi','degradation')

tiledlayout(1,3)
nexttile
fplot(f(1,1,1,1),[0 1])
nexttile
fplot(f(1,1,2,2),[0 1])
nexttile
fplot(f(1,2,1,2),[0 1])

% tiledlayout(1,3)
% nexttile
% fplot(degradation.fun(1,1,1,1),[0 1])
% title('C11')
% nexttile
% fplot(degradation.fun(1,1,2,2),[0 1])
% title('C12')
% nexttile
% fplot(degradation.fun(1,2,1,2),[0 1])
% title('C33')
