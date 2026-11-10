clc,clear,close all
%% Material info
syms x y
alphaAT1 = x;
alphaCOV1 = sqrt(x);
alphaCOV2 = x.^2;

load('HorizontalCrack.mat')
gHC = degradation.fun{2,2,2,2}(x);
load('HCsmlCOV.mat');
CHC_COV = squeeze(Cmat(2,2,2,2,:)); dCHC_COV = squeeze(dCmat(2,2,2,2,:));
load('HCsml.mat')
CHC = squeeze(Cmat(2,2,2,2,:)); dCHC = squeeze(dCmat(2,2,2,2,:));

aPrimeAT1 = double(subs(-diff(alphaAT1),phi));
aPrimeCOV1 = double(subs(-diff(alphaCOV1),phi));
aPrimeCOV2 = double(subs(-diff(alphaCOV2),phi));

gHC_COV = subs(gHC,x,sqrt(y));
aAT1_COV = subs(alphaAT1,x,sqrt(y));

figure(1)
hold on
fplot(gHC,[0 1],'Color','b')
fplot(gHC_COV,[0 1],'LineStyle','-.','Color','b')
plot(phi,CHC,'Color','r')
plot(phi,CHC_COV,'Color','r','LineStyle','-.')
plot(phi.^2,CHC,'Color','g','LineStyle','-.')
plot(sqrt(phi),CHC,'Color','k','LineStyle','-.')
legend({'HC (handle)','HC (handle) $\phi = sqrt(\varphi)$', ...
        'HC (spline)','HC (spline) $\phi = sqrt(\varphi)$',...
        'HC (spline) $f(\phi^2)$','HC (spline) $f(sqrt(\phi))$'},'interpreter','latex')


epsHC  = sqrt(-diff(alphaAT1)/diff(gHC));
sigHC  = sqrt((-diff(alphaAT1)/diff(gHC))*gHC.^2);

epsHC_COV  = sqrt(-diff(aAT1_COV)/diff(gHC_COV));
sigHC_COV  = sqrt((-diff(aAT1_COV)/diff(gHC_COV))*gHC_COV.^2);

epsHC_slm     = sqrt(aPrimeAT1'./dCHC);
epsHC_COV1_slm = sqrt(aPrimeCOV1'./dCHC_COV);
%epsHC_COV2_slm = sqrt(aPrimeCOV2'./dCHC_COV);

sigHC_slm    = sqrt((aPrimeAT1'./dCHC).*CHC.^2);
sigHCCOV_slm = sqrt((aPrimeCOV1'./dCHC_COV).*CHC_COV.^2);


% Plots
figure(4)
hold on
fplot(epsHC,[0 1],'Color','b')
fplot(epsHC_COV,[0 1],'LineStyle','-.','Color','b')
plot(phi,epsHC_slm,'Color','r')
plot(phi,epsHC_COV1_slm,'Color','r','LineStyle','-.')
plot(sqrt(phi),epsHC_slm,'Color','k','LineStyle','-.')
legend({'HC (handle)','HC (handle) $\phi = sqrt(\varphi)$', ...
        'HC (spline)','HC (spline) $\phi = sqrt(\varphi)$',...
        'HC (spline) $f(sqrt(\phi))$'},'interpreter','latex')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')


figure(6)
hold on
%fplot(sigHC,[0 1])
plot(phi,sigHC_slm);
plot(phi,sigHCCOV_slm)
legend('HorizontalCrack (handle)','HorizontalCrack (Spline)','HorizontalCrack COV (Spline)')
title('Stress vs damage')
xlabel('Damage')
ylabel('Stress')




