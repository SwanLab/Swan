clear; close all; clc;

%% INPUT


outDir = 'C:\Users\david\Desktop\PlotsResultats\4.MMB';

if ~exist(outDir,'dir')
mkdir(outDir);
end


%% PLOT SETTINGS

figPos    = [100 100 1200 900];
fontSize  = 22;
lineWidth = 2;

%% MATERIAL PARAMETERS

E  = 120e3;
nu = 0.3;

% E  = E/(1-nu^2);
% nu = nu/(1-nu);

GIc   = 0.260;
GIIc  = 1.002;
eta   = 2;

b     = 1;
h     = 1.55;

a0    = 35;
aMax  = 120;

L     = 75;
c     = 63.18;

%% LEFM SOLUTION

Gc = GIc + (GIIc-GIc)*0.5^eta;

G = E/(2*(1+nu));

Gamma = 1.18*E/G;
chi   = sqrt(E/(11*G)*(3 - 2*(Gamma/(1+Gamma))^2));

a = linspace(a0,aMax,1000);

aModeI  = a + chi*h;
aModeII = a + 0.42*chi*h;

a0ModeI  = a0 + chi*h;
a0ModeII = a0 + 0.42*chi*h;

EI = E*b*h^3/12;

denF = 64*b*L^2*EI;
denU = 96*L^2*EI;

FProp = sqrt(Gc*denF ./ ...
(4*(3*c-L)^2.*aModeI.^2 + 3*(c+L)^2.*aModeII.^2));

uProp = FProp .* ...
(4*(3*c-L)^2.*aModeI.^3 + ...
(c+L)^2.*(3*aModeII.^3 + 2*L^3)) ./ denU;

Fmmb = sqrt(Gc*denF / ...
(4*(3*c-L)^2*a0ModeI^2 + ...
3*(c+L)^2*a0ModeII^2));

Cmmb = (4*(3*c-L)^2*a0ModeI^3 + ...
(c+L)^2*(3*a0ModeII^3 + 2*L^3)) / denU;

ummb = Fmmb*Cmmb;

uElastic = linspace(0,ummb,1000);
FElastic = uElastic/Cmmb;

%% LOAD NUMERICAL RESULTS

% S = load(fileName);
% res15_30 = load('res15_30.mat');
% 
% R{1} = res15_30.res15_30;
% R{2} = S.results(1).output;

R = load("MMBThickLever.mat");
R = R.Res;
labels = {'tau0Shear = 30 MPa',...
          'tau0Shear = 60 MPa'};

markers = {'o','s'};

%% FORCE-DISPLACEMENT PLOT

fig = figure('Position',figPos);

hold on
grid on
box on

plot(uElastic,FElastic,...
    'r',...
    'LineWidth',lineWidth,...
    'DisplayName','LEFM elastic');

plot(uProp,FProp,...
    'r--',...
    'LineWidth',lineWidth,...
    'DisplayName','LEFM propagation');

for i = 1:2

    output = R;
    F = output(:,1);
    U = output(:,2);

    scatter(U,0.5*F,...
        25,...
        markers{i},...
        'filled',...
        'DisplayName',labels{i});

end

xlabel('Lever displacement [mm]','FontSize',fontSize)
ylabel('Load [kN]','FontSize',fontSize)

title('MMB benchmark, $\tau_n^0 = 30$ MPa',...
    'Interpreter','latex',...
    'FontSize',fontSize)
legend('Location','best','FontSize',fontSize-2)

set(gca,'FontSize',fontSize)

xlim([0 10])

exportgraphics(fig,...
    fullfile(outDir,'MMB30Normal.png'),...
    'Resolution',300);