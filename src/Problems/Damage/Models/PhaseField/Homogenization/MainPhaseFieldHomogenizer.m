clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 100;
s.holeType   = 'Hexagon';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = 'Perimeter';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


%% DERIVATIVE
load('SquarePerimeter.mat')
[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('SquarePerimeter','mat','phi','degradation')

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
