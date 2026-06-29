clear;
close all
    
files = {"ENFTauShear30Check.mat","ENFTauShear60.mat","ENFTauShear120.mat"};
labels = {'\tau_t^0 = 30 MPa','\tau_t^0 = 60 MPa','\tau_t^0 = 120 MPa'};
outDir = 'C:\Users\david\Desktop\PlotsResultats\3.ENF';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

figPos = [100 100 1200 900];
fontSize = 22;
lineWidth = 2;

E  = 120e9;
nu = 0.3;
GIIc = 1002;
b  = 1;
h  = 1.55e-3;
a0 = 35e-3;
aMax = 120e-3;
L = 75e-3;

nuStress = nu/(1+nu); 
EStress  = E*(1-nuStress^2);    

E = EStress;
nu = nuStress;

a = linspace(a0,aMax,1000);

% FProp = sqrt((16*b^2*h^3*E*GIIc)./(9*a.^2));
% uProp = FProp.*(3*a.^3 + 2*L.^3)./(8*b*h.^3*E);
% 
% Fenf = sqrt((16*b^2*h^3*E*GIIc)/(9*a0^2));
% Cenf = (3*a0^3 + 2*L^3)/(8*b*h^3*E);
% uenf = Cenf*Fenf;
% 
% uElastic = linspace(0,uenf,1000);
% FElastic = uElastic/Cenf;

G = E/(2*(1+nu));

Gamma = 1.18*E/G;
chi   = sqrt(E/(11*G)*(3 - 2*(Gamma/(1+Gamma))^2));

aEff  = a  + 0.42*chi*h;
a0Eff = a0 + 0.42*chi*h;

FProp = sqrt((16*b^2*h^3*E*GIIc)./(9*aEff.^2));
uProp = FProp.*(3*aEff.^3 + 2*L.^3)./(8*b*h^3*E);

Fenf = sqrt((16*b^2*h^3*E*GIIc)/(9*a0Eff^2));
Cenf = (3*a0Eff^3 + 2*L^3)/(8*b*h^3*E);
uenf = Cenf*Fenf;

uElastic = linspace(0,uenf,1000);
FElastic = uElastic/Cenf;



fig = figure('Position',figPos);

hold on
grid on

plot(uElastic*1e3,FElastic,'r','LineWidth',2,'DisplayName','LEFM elastic')
plot(uProp*1e3,FProp,'r--','LineWidth',2,'DisplayName','LEFM propagation')

markers = {'o','s','d'};

for i = 1:numel(files)
    S = load(files{i});
    fn = fieldnames(S);
    data = S.(fn{1});

    U = data(:,2)*1e3;
    F = data(:,1);

scatter(U,F,30,markers{i},'filled','DisplayName',labels{i})
end

xlabel('Midpoint displacement [mm]')
ylabel('Load [N]')
legend('Location','best')
xlim([0,9])
title("ENF Benchmark. Traction-Displacement Graph")
set(gca,'FontSize',fontSize)



exportgraphics(fig,fullfile(outDir,'ENF.png'),'Resolution',300);