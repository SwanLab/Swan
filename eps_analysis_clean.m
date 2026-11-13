clc,clear,close all
%% Material info
syms x
% Material properties
E = 210; nu = 0.3; 
C22 = E/((1+nu)*(1-nu));
l0 = 0.1; Gc = 5e-3;

% Degradation functions
gAT     = (1-x)^2; 
gAT1COV = (1-x^2)^2;
gAT2COV = (1-sqrt(x))^2;

gLinearSB = 1 - x; % Linear degradation causing limit snap-back
gQuadSB   = 1 - x^2; % Quadratic degradation causing limit snap-back

% Dissipation functions
alphaAT1 = x;          cwAT1 = 8/3; wAT1  = Gc/(l0*cwAT1);
alphaAT2 = x.^2;       cwAT2 = 2;   wAT2  = Gc/(l0*cwAT2);
alphaAT3 = x.^3;       cwAT3 = 8/5; wAT3  = Gc/(l0*cwAT3);
alphaQuad = 2*x - x^2; cwQuad = pi; wQuad = Gc/(l0*cwQuad);
idx = 2; alphaCOV  = x.^(1/idx);
%[alphaHybridCubic, wHybridCubic] = createHybridFunction('Cubic','false',0.3,[-2 2])
[alphaHybridCubic, cwHybridCubic] = createHybridFunction('Cubic','true',0.3,[-0.5]);
%[alphaHybridQuad, wHybridQuad] = createHybridFunction('Quad','false',0.1,[1.8])
[alphaHybridQuad, cwHybridQuad] = createHybridFunction('Quad','true',0.1,[]);
wHybridCubic = Gc/(l0*cwHybridCubic);
wHybridQuad = Gc/(l0*cwHybridQuad);

% Homogenized models (Y-direction -> C(2,2,2,2) ) 
load('SquareArea.mat')
gSqr = degradation.fun{2,2,2,2}(x); dgSqr = degradation.dfun{2,2,2,2}(x);
load('CircleArea.mat') 
gCir = degradation.fun{2,2,2,2}(x); dgCir = degradation.dfun{2,2,2,2}(x);
load('HorizontalCrack.mat')
gHC = degradation.fun{2,2,2,2}(x); phiVHC = phi;
load('CircleCrack.mat')
gCC = degradation.fun{2,2,2,2}(x); phiVCC = phi;

load('HCsml.mat')
CHC = E*squeeze(Cmat(2,2,2,2,:)); dCHC = E*squeeze(dCmat(2,2,2,2,:));
load('HCsmlCOV.mat'); %(phi^2 -> idx = 2)
CHCCOV = E.*squeeze(Cmat(2,2,2,2,:)); dCHCCOV = E.*squeeze(dCmat(2,2,2,2,:));
load('CircleCrack_sml.mat')
CCC = E.*squeeze(Cmat(2,2,2,2,:)); dCCC = E.*squeeze(dCmat(2,2,2,2,:));


%% Strain/stress as a function of damage (all handle)
epsAT1 = sqrt(2*wAT1/C22)*sqrt(-diff(alphaAT1)/diff(gAT));
epsAT2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gAT));
epsAT1COV = sqrt(2*wAT1/C22)*sqrt(-diff(alphaAT2)/diff(gAT1COV));
epsAT2COV = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT1)/diff(gAT2COV));
epsHC  = sqrt(2*wAT1/C22)*sqrt(-diff(alphaAT1)/diff(gHC));  % AT1 dissipation
epsSqr = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gSqr)); % AT2 dissipation 
epsCir = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gCir)); % AT2 dissipation

sigAT1 = sqrt(2*wAT1*C22)*sqrt((-diff(alphaAT1)/diff(gAT))*gAT.^2);
sigAT2 = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gAT))*gAT.^2);
sigAT1COV = sqrt(2*wAT1*C22)*sqrt((-diff(alphaAT2)/diff(gAT1COV))*gAT1COV.^2);
sigAT2COV = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT1)/diff(gAT2COV))*gAT2COV.^2);
sigHC  = sqrt(2*wAT1*C22)*sqrt((-diff(alphaAT1)/diff(gHC))*gHC.^2);
sigSqr = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gSqr))*gSqr.^2);
sigCir = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gCir))*gCir.^2);

% Plots
figure(1)
hold on
fplot(epsAT1,[0 1])
fplot(epsAT2,[0 1])
fplot(epsAT1COV,[0 1])
fplot(epsAT2COV,[0 1])
legend('AT1','AT2','AT1COV','AT2COV')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')

figure(2)
hold on
fplot(epsHC,[0 1])
fplot(epsSqr,[0 1])
fplot(epsCir,[0 1])
legend('HorizontalCrack','Square','Circle')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')

figure(3)
hold on
fplot(sigAT1,[0 1])
fplot(sigAT2,[0 1])
fplot(sigAT2COV,[0 1])
fplot(sigHC,[0 1])
fplot(sigSqr,[0 1])
fplot(sigCir,[0 1])
legend('AT1','AT2','AT2COV','HorizontalCrack','Square','Circle')
title('Stress vs damage')
xlabel('Damage')
ylabel('Stress')


%% Comparing Horizontal Crack (handle vs spline)
aPrime    = double(subs(-diff(alphaAT1),phiVHC));
aPrimeCOV = double(subs(-diff(alphaCOV),phiVHC));

epsHC_slm = sqrt(2*wAT1).*sqrt(aPrime'./dCHC);
epsHCCOV_slm = sqrt(2*wAT1).*sqrt(aPrimeCOV'./dCHCCOV);

sigHC_slm = sqrt(2*wAT1).*sqrt((aPrime'./dCHC).*CHC.^2);
sigHCCOV_slm = sqrt(2*wAT1)*sqrt((aPrimeCOV'./dCHCCOV).*CHCCOV.^2);

% Plots
figure(4)
hold on
fplot(epsHC,[0 1])
plot(phiVHC,epsHC_slm)
plot(phiVHC,epsHCCOV_slm)
legend('HorizontalCrack (handle)','HorizontalCrack (Spline)','HorizontalCrack COV (Spline)')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')

figure(5)
hold on
fplot(sigHC,[0 1])
plot(phiVHC,sigHC_slm);
plot(phiVHC,sigHCCOV_slm)
legend('HorizontalCrack (handle)','HorizontalCrack (Spline)','HorizontalCrack COV (Spline)')
title('Stress vs damage')
xlabel('Damage')
ylabel('Stress')


%% Full brittle response (snap-back limit models)
epsBrittle1 = sqrt(2*wAT1/C22)*sqrt(-diff(alphaAT1)/diff(gLinearSB));
epsBrittle2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gQuadSB));

sigBrittle1 = sqrt(2*wAT1*C22)*sqrt((-diff(alphaAT1)/diff(gLinearSB))*gLinearSB.^2);
sigBrittle2 = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gQuadSB))*gQuadSB.^2);

% figure(6)
% hold on
% fplot(epsBrittle1,[0 1])
% fplot(epsBrittle2,[0 1],'LineStyle','--')
% legend({'$g(\phi) = 1-\phi^2; \omega(\phi) = \phi^2$','$g(\phi) = 1-\phi; \omega(\phi) = \phi$',...
%         },'Interpreter','latex')
% title('Strain vs damage')
% xlabel('Damage')
% ylabel('Strain')


% figure(7)
% hold on
% fplot(sigBrittle1,[0 1])
% fplot(sigBrittle2,[0 1],'LineStyle','--')
% legend({'$g(\phi) = 1-\phi^2; \omega(\phi) = \phi^2$','$g(\phi) = 1-\phi; \omega(\phi) = \phi$',...
%         },'Interpreter','latex')
% title('Stress vs damage')
% xlabel('Damage')
% ylabel('Stress')

%% Horizontal crack with different dissipations
aPrimeAT2  = double(subs(-diff(alphaAT2),phiVHC));
aPrimeAT3  = double(subs(-diff(alphaAT3),phiVHC));
aPrimeQuad = double(subs(-diff(alphaQuad),phiVHC));
aPrimeHybridCubic = double(subs(-diff(alphaHybridCubic),phiVHC)); 
aPrimeHybridQuad  = double(subs(-diff(alphaHybridQuad),phiVHC));

epsHC1 = sqrt(2*wAT1/C22)*sqrt(-diff(alphaAT1)/diff(gHC));
epsHC2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gHC));
epsHC3 = sqrt(2*wAT3/C22)*sqrt(-diff(alphaAT3)/diff(gHC));
epsHCQuad = sqrt(2*wQuad/C22)*sqrt(-diff(alphaQuad)/diff(gHC));
epsHCHybridQuad  = sqrt(2*wHybridQuad/C22)*sqrt(-diff(alphaHybridQuad)/diff(gHC));
epsHCHybridCubic = sqrt(2*wHybridCubic/C22)*sqrt(-diff(alphaHybridCubic)/diff(gHC));

epsHC1_slm = sqrt(2*wAT1).*sqrt(aPrime'./dCHC);
epsHC2_slm = sqrt(2*wAT2).*sqrt(aPrimeAT2'./dCHC);
epsHC3_slm = sqrt(2*wAT3).*sqrt(aPrimeAT3'./dCHC);
epsHCquad_slm = sqrt(2*wQuad).*sqrt(aPrimeQuad'./dCHC);
epsHCHybridQuad_slm = sqrt(2*wHybridQuad)*sqrt(aPrimeHybridQuad'./dCHC);
epsHCHybridCubic_slm = sqrt(2*wHybridCubic)*sqrt(aPrimeHybridCubic'./dCHC);

% Plots
figure(8)
hold on
fplot(alphaAT1,[0 1])
fplot(alphaAT2,[0 1])
fplot(alphaAT3,[0 1])
fplot(alphaQuad,[0 1])
fplot(alphaHybridQuad,[0 1])
fplot(alphaHybridCubic,[0 1])
title('Horizontal crack dissipation functions')
xlabel('Damage')
ylabel('Dissipation')
legend('Linear','Quadratic','Cubic','Quadratic (cohesive)', ...
       'Quadratic - Linear','Cubic - Linear')

figure(9)
hold on
fplot(epsHC1,[0 1])
fplot(epsHC2,[0 1])
fplot(epsHC3,[0 1])
fplot(epsHCQuad,[0 1])
fplot(epsHCHybridQuad,[0 1])
fplot(epsHCHybridCubic,[0 1])
legend({'HC linear','HC quadratic','HC cubic','HC quadratic (cohesive)', ...
        'HC hybrid Quad','HC hybrid Cubic'},'Interpreter','latex')
title('Strain vs damage (handle)')
xlabel('Damage')
ylabel('Strain')

figure(10)
hold on
plot(phiVHC,epsHC1_slm)
plot(phiVHC,epsHC2_slm)
plot(phiVHC,epsHC3_slm)
plot(phiVHC,epsHCquad_slm)
plot(phiVHC,epsHCHybridQuad_slm)
plot(phiVHC,epsHCHybridCubic_slm)
legend({'HC linear','HC quadratic','HC cubic','HC quadratic (cohesive)', ...
        'HC hybrid Quad','HC hybrid Cubic'},'Interpreter','latex')
title('Strain vs damage (spline)')
xlabel('Damage')
ylabel('Strain')

%% Circular crack response 
[alphaHybridCubic, cwHybridCubic] = createHybridFunction('Cubic','false',0.05,[5 5]);
[alphaHybridQuad, cwHybridQuad] = createHybridFunction('Quad','false',0.05,[5]);
wHybridCubic = Gc/(l0*cwHybridCubic);
wHybridQuad = Gc/(l0*cwHybridQuad);

aPrime    = double(subs(-diff(alphaAT1),phiVCC));
aPrimeAT2  = double(subs(-diff(alphaAT2),phiVCC));
aPrimeAT3  = double(subs(-diff(alphaAT3),phiVCC));
aPrimeQuad = double(subs(-diff(alphaQuad),phiVCC));
aPrimeHybridCubic = double(subs(-diff(alphaHybridCubic),phiVCC));
aPrimeHybridQuad  = double(subs(-diff(alphaHybridQuad),phiVCC));  

epsCC1 = sqrt(2*wAT1/C22)*sqrt(-diff(alphaAT1)/diff(gCC));
epsCC2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gCC));
epsCC3 = sqrt(2*wAT3/C22)*sqrt(-diff(alphaAT3)/diff(gCC));
epsCCQuad = sqrt(2*wQuad/C22)*sqrt(-diff(alphaQuad)/diff(gCC));
epsCCHybridQuad  = sqrt(2*wHybridQuad/C22)*sqrt(-diff(alphaHybridQuad)/diff(gCC));
epsCCHybridCubic = sqrt(2*wHybridCubic/C22)*sqrt(-diff(alphaHybridCubic)/diff(gCC));

epsCC1_slm = sqrt(2*wAT1).*sqrt(aPrime./dCCC);
epsCC2_slm = sqrt(2*wAT2).*sqrt(aPrimeAT2./dCCC);
epsCC3_slm = sqrt(2*wAT3).*sqrt(aPrimeAT3./dCCC);
epsCCquad_slm = sqrt(2*wQuad).*sqrt(aPrimeQuad./dCCC);
epsCCHybridQuad_slm = sqrt(2*wHybridQuad)*sqrt(aPrimeHybridQuad./dCCC);
epsCCHybridCubic_slm = sqrt(2*wHybridCubic)*sqrt(aPrimeHybridCubic./dCCC);

figure(11)
hold on
fplot(alphaAT1,[0 1])
fplot(alphaAT2,[0 1])
fplot(alphaAT3,[0 1])
fplot(alphaQuad,[0 1])
fplot(alphaHybridQuad,[0 1])
fplot(alphaHybridCubic,[0 1])
title('Circle + crack dissipation functions')
xlabel('Damage')
ylabel('Dissipation')
legend('Linear','Quadratic','Cubic','Quadratic (cohesive)', ...
       'Quadratic - Linear','Cubic - Linear')

figure(12)
hold on
fplot(epsCC1,[0 1])
fplot(epsCC2,[0 1])
fplot(epsCC3,[0 1])
fplot(epsCCQuad,[0 1])
fplot(epsCCHybridQuad,[0 1])
fplot(epsCCHybridCubic,[0 1])
legend({'CC linear','CC quadratic','CC cubic','CC quadratic (cohesive)', ...
        'CC hybrid Quad','CC hybrid Cubic'},'Interpreter','latex')
title('Strain vs damage (handle)')
xlabel('Damage')
ylabel('Strain')

figure(13)
hold on
plot(phiVCC,epsCC1_slm)
plot(phiVCC,epsCC2_slm)
plot(phiVCC,epsCC3_slm)
plot(phiVCC,epsCCquad_slm)
plot(phiVCC,epsCCHybridQuad_slm)
plot(phiVCC,epsCCHybridCubic_slm)
legend({'CC linear','CC quadratic','CC cubic','CC quadratic (cohesive)', ...
        'CC hybrid Quad','CC hybrid Cubic'},'Interpreter','latex')
title('Strain vs damage (spline)')
xlabel('Damage')
ylabel('Strain')











