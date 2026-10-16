close all
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
load('HorizontalCrack.mat')
gHC = degradation.fun{2,2,2,2}(x);
load('SquareArea.mat')
gSqr = degradation.fun{2,2,2,2}(x); dgSqr = degradation.dfun{2,2,2,2}(x);
load('CircleArea.mat') 
gCir = degradation.fun{2,2,2,2}(x); dgCir = degradation.dfun{2,2,2,2}(x);

% Homogenized crack with splines
load('HCsmlCOV.mat');
CmatCOV = 210.*Cmat; dCmatCOV = 210.*dCmat; d2CmatCOV = 210.*d2Cmat;
load('HCsml.mat')
Cmat = 210.*Cmat; dCmat = 210.*dCmat; d2Cmat = 210.*d2Cmat;

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
sigAT2COV = sqrt(2*wAT2*C22)*sqrt((-diff(alpha)/diff(gATCOV))*gATCOV.^2);
sigHC  = sqrt(2*w*C22)*sqrt((-diff(alpha)/diff(gHC))*gHC.^2);
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

% figure(3)
% hold on
% fplot(sigAT1,[0 1])
% fplot(sigAT2,[0 1])
% fplot(sigAT2COV,[0 1])
% fplot(sigHC,[0 1])
% fplot(sigSqr,[0 1])
% fplot(sigCir,[0 1])
% legend('AT1','AT2','AT2COV','HorizontalCrack','Square','Circle')
% title('Stress vs damage')
% xlabel('Damage')
% ylabel('Stress')


%% Strain as a function of damage (points [SLM])



alphaProposed = x.^2;
aPrime = double(subs(-diff(alpha),phi));
aPrimeCOV = double(subs(-diff(alphaProposed),phi));

epsHC = sqrt(2*w).*sqrt(aPrimeCOV'./squeeze(dCmat(2,2,2,2,:)));
sigHomog = sqrt(2*w).*sqrt((aPrimeCOV'./squeeze(dCmat(2,2,2,2,:))).*squeeze(Cmat(2,2,2,2,:)).^2);
sigHomogCOV = sqrt(2*w)*sqrt((aPrimeCOV'./squeeze(dCmatCOV(2,2,2,2,:))).*squeeze(CmatCOV(2,2,2,2,:)).^2);


% Plots
figure(4)
hold on
plot(phi,epsHomog)
plot(phi,epsHomogCOV)
legend('HorizontalCrack (handle)','HorizontalCrack (Spline)')
title('Strain vs damage')
xlabel('Damage')
ylabel('Strain')


figure(6)
hold on
plot(phi,sigHomog);
plot(phi,epsHomogCOV)
legend('HorizontalCrack (handle)','HorizontalCrack (Spline)')
title('Stress vs damage')
xlabel('Damage')
ylabel('Stress')

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



