clc,clear,close all
%% Material info
% AT model 1D info
syms x
gAT = (1-x)^2;
gAT2COV = (1-sqrt(x))^2;
gAT1COV = (1-x^2)^2;
alpha = x; cw = 8/3;
alphaAT2 = x.^2; cwAT2 = 2;

E = 210; nu = 0.3; 
C22 = E/((1+nu)*(1-nu));
l0 = 0.1; Gc = 5e-3;
w = Gc/(l0*cw);
wAT2 = Gc/(l0*cwAT2);

% Homogenized model 1D info (Y-direction)
load('SquareArea.mat')
gSqr = degradation.fun{2,2,2,2}(x); dgSqr = degradation.dfun{2,2,2,2}(x);
load('CircleArea.mat') 
gCir = degradation.fun{2,2,2,2}(x); dgCir = degradation.dfun{2,2,2,2}(x);

load('HorizontalCrack.mat')
gHC = degradation.fun{2,2,2,2}(x); phiV = phi;
load('HCsmlCOV.mat');
CHCCOV = 210.*squeeze(Cmat(2,2,2,2,:)); dCHCCOV = 210.*squeeze(dCmat(2,2,2,2,:));
load('HCsml.mat')
CHC = 210.*squeeze(Cmat(2,2,2,2,:)); dCHC = 210.*squeeze(dCmat(2,2,2,2,:));
load('CircleCrack.mat')
gCC = degradation.fun{2,2,2,2}(x); phiVCC = phi;
load('CircleCrack_sml.mat')
CCC = 210.*squeeze(Cmat(2,2,2,2,:)); dCCC = 210.*squeeze(dCmat(2,2,2,2,:));


%% Strain as a function of damage (handle)
epsAT1 = sqrt(2*w/C22)*sqrt(-diff(alpha)/diff(gAT));
epsAT2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gAT));
epsAT1COV = sqrt(2*w/C22)*sqrt(-diff(alphaAT2)/diff(gAT1COV));
epsAT2COV = sqrt(2*wAT2/C22)*sqrt(-diff(alpha)/diff(gAT2COV));
epsHC  = sqrt(2*w/C22)*sqrt(-diff(alpha)/diff(gHC));
epsSqr = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gSqr));
epsCir = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gCir));

sigAT1 = sqrt(2*w*C22)*sqrt((-diff(alpha)/diff(gAT))*gAT.^2);
sigAT2 = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gAT))*gAT.^2);
sigAT1COV = sqrt(2*w*C22)*sqrt((-diff(alphaAT2)/diff(gAT1COV))*gAT1COV.^2);
sigAT2COV = sqrt(2*wAT2*C22)*sqrt((-diff(alpha)/diff(gAT2COV))*gAT2COV.^2);
sigHC  = sqrt(2*w*C22)*sqrt((-diff(alpha)/diff(gHC))*gHC.^2);
sigSqr = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gSqr))*gSqr.^2);
sigCir = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gCir))*gCir.^2);

% 
% epsAT2new = subs(epsAT2,x,sqrt(y));
% dEPsAT2COV = diff(epsAT2COV);
% dEPsAT2new = diff(epsAT2new);
% dEPsAT2computed = subs(diff(epsAT2),x,sqrt(y))*(0.5/sqrt(y));
% dEPsAT2computed = diff(epsAT2)*(0.5/x);
% 
% figure()
% hold on
% fplot(epsAT2COV,[0 1])
% fplot(epsAT2new,[0 1],'LineStyle','--')
% 
% figure()
% hold on
% fplot(dEPsAT2COV,[0 1])
% fplot(dEPsAT2new,[0 1],'LineStyle','--')
% fplot(dEPsAT2computed,[0 1],'LineStyle',':','Color','k');


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


%% Strain as a function of damage (points [SLM])

% Homogenized crack with splines
idx = 2;
alphaCOV  = x.^(1/idx);
aPrime    = double(subs(-diff(alpha),phiV));
aPrimeCOV = double(subs(-diff(alphaCOV),phiV));

epsHC_slm = sqrt(2*w).*sqrt(aPrime'./dCHC);
epsHCCOV_slm = sqrt(2*w).*sqrt(aPrimeCOV'./dCHCCOV);

sigHC_slm = sqrt(2*w).*sqrt((aPrime'./dCHC).*CHC.^2);
sigHCCOV_slm = sqrt(2*w)*sqrt((aPrimeCOV'./dCHCCOV).*CHCCOV.^2);


% Plots
figure(4)
hold on
%fplot(epsHC,[0 1])
plot(phiV,epsHC_slm)
plot(phiV,epsHCCOV_slm)
legend('HorizontalCrack (handle)','HorizontalCrack (Spline)','HorizontalCrack COV (Spline)')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')


figure(6)
hold on
%fplot(sigHC,[0 1])
plot(phiV,sigHC_slm);
plot(phiV,sigHCCOV_slm)
legend('HorizontalCrack (handle)','HorizontalCrack (Spline)','HorizontalCrack COV (Spline)')
title('Stress vs damage')
xlabel('Damage')
ylabel('Stress')


%% SLM plots
close all
gB1 = 1 - x;
gB2 = 1 - x^2;
alphaAT3 = x.^3; cwAT3 = 8/5;
alphaQuad = 2*x - x^2; cwQuad = pi;
wAT3 = Gc/(l0*cwAT3); wQuad = Gc/(l0*cwQuad);

aPrimeAT2 = double(subs(-diff(alphaAT2),phiV));
aPrimeAT3 = double(subs(-diff(alphaAT3),phiV));
aPrimeQuad = double(subs(-diff(alphaQuad),phiV));


% a = -2; b = 2; x0 = 0.3;
% c = 2*a*x0^3 + (b-3*a)*x0^2 - 2*b*x0 + 1;
% m = 3*a*x0^2 + 2*b*x0 + c;
% alphaHybridCubic = piecewise((x<x0),(a*x^3 + b*x^2 + c*x),(x>=x0),(m*x-m+1));
% aPrimeHybridCubic = double(subs(-diff(alphaHybridCubic),phiV)); 
% wHybridCubic = double(4*int(sqrt(alphaHybridCubic),0,1));
% 
% a = 1.8; x0 = 0.3;
% b = a*x0^2 - 2*a*x0 + 1;
% m = 2*a*x0 + b;
% alphaHybridQuad = piecewise((x<x0),(a*x^2 + b*x),(x>=x0),(m*x-m+1));
% aPrimeHybridQuad = double(subs(-diff(alphaHybridQuad),phiV)); 
% wHybridQuad = double(4*int(sqrt(alphaHybridQuad),0,1));

a = -0.5; x0 = 0.3;
b = (2*a*x0^3 - 3*a*x0^2 + 1)/(-x0^2 + 2*x0); 
m = 3*a*x0^2 + 2*b*x0;
alphaHybridCubic = piecewise((x<x0),(a*x^3 + b*x^2),(x>=x0),(m*x-m+1));
aPrimeHybridCubic = double(subs(-diff(alphaHybridCubic),phiV)); 
wHybridCubic = double(4*int(sqrt(alphaHybridCubic),0,1));

x0 = 0.1; a = -1/(x0^2 - 2*x0);
m = 2*a*x0;
alphaHybridQuad = piecewise((x<x0),(a*x^2),(x>=x0),(m*x-m+1));
aPrimeHybridQuad = double(subs(-diff(alphaHybridQuad),phiV)); 
wHybridQuad = double(4*int(sqrt(alphaHybridQuad),0,1));

epsBrittle1 = sqrt(2*w/C22)*sqrt(-diff(alpha)/diff(gB1));
epsBrittle2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gB2));
epsHC1_slm = sqrt(2*w).*sqrt(aPrime'./dCHC);
epsHC2_slm = sqrt(2*wAT2).*sqrt(aPrimeAT2'./dCHC);
epsHC3_slm = sqrt(2*wAT3).*sqrt(aPrimeAT3'./dCHC);
epsHCquad_slm = sqrt(2*wQuad).*sqrt(aPrimeQuad'./dCHC);
epsHCHybridCubic_slm = sqrt(2*wHybridCubic/C22)*sqrt(aPrimeHybridCubic'./dCHC);
epsHCHybridQuad_slm = sqrt(2*wHybridQuad/C22)*sqrt(aPrimeHybridQuad'./dCHC);

epsHC1 = sqrt(2*w/C22)*sqrt(-diff(alpha)/diff(gHC));
epsHC2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gHC));
epsHC3 = sqrt(2*wAT3/C22)*sqrt(-diff(alphaAT3)/diff(gHC));
epsHCHybrid = sqrt(2*wHybridCubic/C22)*sqrt(-diff(alphaHybridCubic)/diff(gHC));

sigBrittle1 = sqrt(2*w*C22)*sqrt((-diff(alpha)/diff(gB1))*gB1.^2);
sigBrittle2 = sqrt(2*wAT2*C22)*sqrt((-diff(alphaAT2)/diff(gB2))*gB2.^2);
sigHC1 = sqrt(2*w).*sqrt((aPrime'./dCHC).*CHC.^2);
sigHC2 = sqrt(2*wAT2).*sqrt((aPrimeAT2'./dCHC).*CHC.^2);
sigHC3 = sqrt(2*wAT3).*sqrt((aPrimeAT3'./dCHC).*CHC.^2);


% Plots
figure(100)
hold on
fplot(alpha,[0 1])
fplot(alphaAT2,[0 1])
fplot(alphaAT3,[0 1])
fplot(alphaQuad,[0 1])
fplot(alphaHybridQuad,[0 1])
fplot(alphaHybridCubic,[0 1])
xlabel('Damage')
ylabel('Dissipation')
legend('Linear','Quadratic','Cubic','Quadratic (cohesive)','Quadratic - Linear','Cubic - Linear')

% figure(7)
% hold on
% fplot(epsBrittle1,[0 1])
% fplot(epsBrittle2,[0 1],'LineStyle','--')
% legend({'$g(\phi) = 1-\phi^2; \omega(\phi) = \phi^2$','$g(\phi) = 1-\phi; \omega(\phi) = \phi$',...
%         },'Interpreter','latex')
% title('Strain vs damage')
% xlabel('Damage')
% ylabel('Strain')


% figure(8)
% hold on
% fplot(sigBrittle1,[0 1])
% fplot(sigBrittle2,[0 1],'LineStyle','--')
% legend({'$g(\phi) = 1-\phi^2; \omega(\phi) = \phi^2$','$g(\phi) = 1-\phi; \omega(\phi) = \phi$',...
%         },'Interpreter','latex')
% title('Stress vs damage')
% xlabel('Damage')
% ylabel('Stress')

figure(9)
hold on
plot(phiV,epsHC1_slm)
plot(phiV,epsHC2_slm)
plot(phiV,epsHC3_slm)
plot(phiV,epsHCquad_slm)
plot(phiV,epsHCHybridQuad_slm)
plot(phiV,epsHCHybridCubic_slm)
legend({'HC linear','HC quadratic','HC cubic','HC quadratic (cohesive)','HC hybrid Quad','HC hybrid Cubic'},'Interpreter','latex')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')

% figure(10)
% hold on
% plot(phiV,sigHC1)
% plot(phiV,sigHC2)
% plot(phiV,sigHC3)
% legend({'HC linear','HC quadratic','HC cubic'},'Interpreter','latex')
% title('Stress vs damage')
% xlabel('Damage')
% ylabel('Stress')

% figure(11)
% hold on
% fplot(epsHC1,[0 1])
% fplot(epsHC2,[0 1])
% fplot(epsHC3,[0 1])
% fplot(epsHCHybrid,[0 1])
% legend({'HC linear','HC quadratic','HC cubic','HC hybrid'},'Interpreter','latex')
% title('Strain vs damage')
% xlabel('Damage')
% ylabel('Strain')

%% Circular crack vs Horizontal crack
epsHC  = sqrt(2*w/C22)*sqrt(-diff(alpha)/diff(gHC));
epsCC  = sqrt(2*w/C22)*sqrt(-diff(alpha)/diff(gCC));
epsCCAT2 = sqrt(2*wAT2/C22)*sqrt(-diff(alphaAT2)/diff(gCC));
sigHC  = sqrt(2*w*C22)*sqrt((-diff(alpha)/diff(gHC))*gHC.^2);
sigCC  = sqrt(2*w*C22)*sqrt((-diff(alpha)/diff(gCC))*gCC.^2);

a = 1; x0 = 0.1;
m = 2*a*x0;
alphaHybridQuad = piecewise((x<x0),(a*x^2),(x>=x0),(m*x-m+1));
aPrimeHybridQuad = double(subs(-diff(alphaHybridQuad),phiVCC)); 
wHybridQuad = double(4*int(sqrt(alphaHybridQuad),0,1));


aPrimeQuad = double(subs(-diff(alphaQuad),phiVCC)); 
wQuad = double(4*int(sqrt(alphaQuad),0,1));

gCCv = double(subs(gCC,phiVCC));
epsCCquad = sqrt(2*wQuad/C22)*sqrt(aPrimeQuad'./gCCv);


figure(12)
hold on
fplot(epsHC,[0 1])
fplot(epsCC,[0 1])
fplot(epsCCAT2,[0 1])
legend('HorizontalCrack','Circle + crack','CC AT2')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')

% figure(13)
% hold on
% fplot(sigHC,[0 1])
% fplot(sigCC,[0 1])
% legend('HorizontalCrack','Circle + crack')
% title('Stress vs damage')
% xlabel('Damage')
% ylabel('Stress')














%% Computation 2D (Isotropic materials)
% syms x 
% [dataHomog]  = load('HorizontalCrack.mat');
% C11Homog = dataHomog.degradation.fun{1,1,1,1}(x);
% C12Homog = dataHomog.degradation.fun{2,2,1,1}(x);
% C33Homog = dataHomog.degradation.fun{1,2,1,2}(x);
% 
% bulkHomog   = (C11Homog-C33Homog);
% shearHomog  = C33Homog;
% 
% 
% k = bulkHomog;
% m = shearHomog;
% CPrime1DHomog = ( (diff(k)+diff(m))*(1+((k-m)/(k+m))^2) - 2*(diff(k)-diff(m))*((k-m)/(k+m)));
% 
% 
% epsHomog = sqrt(-diff(wAT1)/CPrime1DHomog);
% figure(102)
% a = fplot(epsHomog,[0 1]);
% aXData = a.XData; aYData = a.YData;
% plot(aYData,aXData)
% legend('Homog')
% title('Damage-displacement')
% 
% figure(103)
% hold on
% fplot(CPrime1DHomog,[0 1])
% legend('Homog')
% title('Degradation derivative equivalent to 1D')



