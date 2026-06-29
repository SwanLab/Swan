clear;

files  = {'DCBTauNormal15.mat','DCBTauNormal30.mat','DCBTauNormal60Check.mat'};
labels = {'\tau_n^0 = 15 MPa','\tau_n^0 = 30 MPa','\tau_n^0 = 60 MPa'};

outDir = 'C:\Users\david\Desktop\PlotsResultats\2.DCB';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

figPos = [100 100 1200 900];
fontSize = 22;
lineWidth = 2;

figPos = [100 100 1200 900];
fontSize = 22;
lineWidth = 2;

E  = 120e9;
nu = 0.3;
Gc = 260;
b  = 1;
h  = 1.55e-3;
a0 = 35e-3;
aMax = 80e-3;

nuStress = nu/(1+nu); 
EStress  = E*(1-nuStress^2);    

E = EStress;
nu = nuStress;

G13   = E/(2*(1+nu));
Gamma = 1.18*sqrt(E*E)/G13;
chi   = sqrt((E/(11*G13))*(3-2*(Gamma/(1+Gamma))^2));

Fdcb = sqrt(E*b^2*h^3*Gc/(12*(a0+chi*h)^2));
Cdcb = 8*(a0+chi*h)^3/(E*b*h^3);

uElastic = linspace(0,3e-3,1000);
FElastic = uElastic/Cdcb;

a = linspace(a0*0.9,aMax,1000);
FProp = sqrt(E*b^2*h^3*Gc./(12*(a+chi*h).^2));
CProp = 8*(a+chi*h).^3./(E*b*h^3);
uProp = CProp.*FProp;

fig = figure('Position',figPos);
hold on
grid on

plot(uElastic*1e3,FElastic,'r','LineWidth',lineWidth,'DisplayName','LEFM elastic')
plot(uProp*1e3,FProp,'r--','LineWidth',lineWidth,'DisplayName','LEFM propagation')

markers = {'o','s','d'};

for i = 1:numel(files)
    S = load(files{i});
    fn = fieldnames(S);
    data = S.(fn{1});

    U = 2*data(:,2)*1e3;
    F = data(:,1);

    scatter(U,F,30,markers{i},'filled','DisplayName',labels{i})
end

xlabel('Opening displacement [mm]','FontSize',fontSize)
ylabel('Load [N]','FontSize',fontSize)
title('DCB Benchmark. Traction displacement graph','FontSize',fontSize)
legend('Location','best','FontSize',fontSize-2)
set(gca,'FontSize',fontSize)

xlim([0,6])


exportgraphics(fig,fullfile(outDir,'DCBPlaneStress.png'),'Resolution',300);