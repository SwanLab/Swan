close all
E  = 1;
%nu = 0.3;
k  = @(nu) E./(2.*(1-nu));
mu = @(nu) E./(2.*(1+nu));
C11 = @(nu) E/((1+nu)*(1-nu));

k0     = 1e-10;
k1     = @(nu) k(nu);
mu0    = 1e-10;
mu1    = @(nu) mu(nu);
etak0  = mu0;
etak1  = @(nu) mu1(nu);
etamu0 = (k0.*mu0)./(2.*mu0+k0);
etamu1 = @(nu) (k1(nu).*mu1(nu))./(2.*mu1(nu)+k1(nu));

kUB  = @(phi,nu) (1/k1(nu)).*(k0.*(phi) + k1(nu).*(1-phi) - ((1-phi).*phi.*(k1(nu)-k0).^2)./(k0.*(1-phi) + k1(nu).*phi + etak1(nu)));
muUB = @(phi,nu) (1/mu1(nu)).*(mu0.*(phi) + mu1(nu).*(1-phi) - ((1-phi).*phi.*(mu1(nu)-mu0).^2)./(mu0.*(1-phi) + mu1(nu).*phi + etamu1(nu)));
kLB  = @(phi,nu) (1/k1(nu)).*(k0.*(phi) + k1(nu).*(1-phi) - ((1-phi).*phi.*(k1(nu)-k0).^2)./(k0.*(1-phi) + k1(nu).*phi + etak0));
muLB = @(phi,nu) (1/mu1(nu)).*(mu0.*(phi) + mu1(nu).*(1-phi) - ((1-phi).*phi.*(mu1(nu)-mu0).^2)./(mu0.*(1-phi) + mu1(nu).*phi + etamu0));

%% AT1 and AT2 functions
gAT1 = @(phi) (1-phi)^2;
gAT2 = @(phi) (1-sqrt(phi))^2;

%% Figure 1: AT1/AT2 vs H-S bounds
cmp = orderedcolors("gem");
cmpGrad = gray(10);
figure(1)
t = tiledlayout(1,2);
nexttile
hold on
grid minor
fplot(@(phi) kUB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','--','LineWidth',3)
ylabel('$\kappa(\phi)/\kappa_0$ [-]','Interpreter','latex');
xlabel({"$\phi$ [-]";"(a)"},'Interpreter','latex');
nexttile
hold on
grid minor

fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p3 = fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p1 = fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3);
p2 = fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','--','LineWidth',3);
ylabel('$\mu(\phi)/\mu_0$ [-]','Interpreter','latex');
xlabel({"$\phi$ [-]";"(b)"},'Interpreter','latex');

fontsize(gcf,40,'points')
lgd = legend([p1,p2,p3],'AT1','AT2','H-S bounds');
lgd.FontSize = 30;

%% Figure 2: AT vs H-S bounds for different Possion

cmpGrad = gray(10)
figure(2)
t = tiledlayout(1,2);
nexttile
hold on
grid minor
fplot(@(phi) kUB(phi,-0.9),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,-0.5),[0 1],'Color',cmpGrad(2,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0),   [0 1],'Color',cmpGrad(3,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0.3),[0 1],'Color',cmpGrad(4,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0.5),[0 1],'Color',cmpGrad(5,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(6,:),'LineStyle','-','LineWidth',3);
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','--','LineWidth',3)
ylabel('$\kappa(\phi)/\kappa_0$ [-]','Interpreter','latex');
xlabel({'$\phi$ [-]';"(a)"},'Interpreter','latex');

nexttile
hold on
grid minor
p3 = fplot(@(phi) muUB(phi,-0.9),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,-0.5),[0 1],'Color',cmpGrad(2,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,0),   [0 1],'Color',cmpGrad(3,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmpGrad(4,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,0.5),[0 1],'Color',cmpGrad(5,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p1 = fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3);
p2 = fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','--','LineWidth',3);
ylabel('$\mu(\phi)/\mu_0$ [-]','Interpreter','latex');
xlabel({"$\phi$ [-]";"(b)"},'Interpreter','latex');

fontsize(gcf,40,'points')
lgd = legend([p1,p2,p3],'AT1','AT2','H-S bounds');
lgd.FontSize = 30;

%% Alessi and Wu functions
gRational = @(phi,gamma) (1-phi)/((1-phi)+gamma*phi);


%% Figure 3_ Rational vs H-S bounds
cmpGamma = autumn(12);

figure(4)
t = tiledlayout(1,2);
nexttile
hold on
grid minor
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) gRational(phi,1),[0 1],'Color',cmpGamma(6,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) gRational(phi,3),[0 1],'Color',cmpGamma(7,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) gRational(phi,5),[0 1],'Color',cmpGamma(8,:),'LineStyle','-','LineWidth',3)
fplot(@(phi) gRational(phi,10),[0 1],'Color',cmpGamma(9,:),'LineStyle','-','LineWidth',3)
fplot(@(phi) gRational(phi,50),[0 1],'Color',cmpGamma(10,:),'LineStyle','-','LineWidth',3);
ylabel('$\kappa(\phi)/\kappa_0$ [-]','Interpreter','latex');
xlabel({'$\phi$ [-]';"(a)"},'Interpreter','latex');

nexttile
hold on
grid minor
p2 = fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) gRational(phi,1),[0 1],'Color',cmpGamma(6,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) gRational(phi,3),[0 1],'Color',cmpGamma(7,:),'LineStyle','-','LineWidth',3)
fplot(@(phi) gRational(phi,5),[0 1],'Color',cmpGamma(8,:),'LineStyle','-','LineWidth',3)
fplot(@(phi) gRational(phi,10),[0 1],'Color',cmpGamma(9,:),'LineStyle','-','LineWidth',3)
p1 = fplot(@(phi) gRational(phi,50),[0 1],'Color',cmpGamma(10,:),'LineStyle','-','LineWidth',3);
ylabel('$\mu(\phi)/\mu_0$ [-]','Interpreter','latex');
xlabel({"$\phi$ [-]";"(b)"},'Interpreter','latex');

fontsize(gcf,40,'points')
lgd = legend([p1,p2],{'Rational','H-S bounds'},'Interpreter','latex');
lgd.FontSize = 30;

%% Figure 4 Constitutive tensor homogenized
[dataHexa]  = load('HexagonBenchmark03.mat');
C11hexaMat = squeeze(dataHexa.mat(1,1,1,1,:));
C12hexaMat = squeeze(dataHexa.mat(2,2,1,1,:));
C33hexaMat = squeeze(dataHexa.mat(1,2,1,2,:));
[dataHoney] = load('HoneycombBenchmark03.mat');
C11honeyMat = squeeze(dataHoney.mat(1,1,1,1,:));
C12honeyMat = squeeze(dataHoney.mat(2,2,1,1,:));
C33honeyMat = squeeze(dataHoney.mat(1,2,1,2,:));

figure(4)
t = tiledlayout(1,3);
nexttile
hold on
grid minor
fplot(dataHexa.degradation.fun{1,1,1,1},[0 1],'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
fplot(dataHoney.degradation.fun{1,1,1,1},[0 1],'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel(char(8450)+"11 [GPa]");
xlabel({"$\phi$ [-]";"(a)"},'Interpreter','latex');
ylim([0,inf])

nexttile
hold on
grid minor
fplot(dataHexa.degradation.fun{2,2,1,1},[0 1],'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
fplot(dataHoney.degradation.fun{2,2,1,1},[0 1],'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel(char(8450)+"12 [GPa]");
xlabel({"$\phi$ [-]";"(b)"},'Interpreter','latex');
ylim([0,inf])
fontsize(gcf,30,'points')
nexttile
hold on
grid minor
fplot(dataHexa.degradation.fun{1,2,1,2},[0 1],'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
fplot(dataHoney.degradation.fun{1,2,1,2},[0 1],'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel(char(8450)+"33 [GPa]");
ylim([0,inf])
xlabel({"$\phi$ [-]";"(c)"},'Interpreter','latex');

fontsize(gcf,40,'points')
lgd = legend('Hexagon','Reinforced hexagon','Orientation','horizontal');
lgd.FontSize = 30;

%% Figure 5 Zener Ratio
C11hexa = dataHexa.degradation.fun{1,1,1,1};
C12hexa = dataHexa.degradation.fun{2,2,1,1};
C33hexa = dataHexa.degradation.fun{1,2,1,2};
C11honey = dataHoney.degradation.fun{1,1,1,1};
C12honey = dataHoney.degradation.fun{2,2,1,1};
C33honey = dataHoney.degradation.fun{1,2,1,2};
ZenerRatioHexa = @(phi) 2*C33hexa(phi)./(C11hexa(phi)-C12hexa(phi));
ZenerRatioHoney = @(phi) 2*C33honey(phi)./(C11honey(phi)-C12honey(phi));

figure(6)
hold on
grid minor
fplot(ZenerRatioHexa,[0 1],'Color',cmp(4,:),'LineStyle','-','LineWidth',3);
fplot(ZenerRatioHoney,[0 1],'Color',cmp(4,:),'LineStyle','--','LineWidth',3);
ylabel("Zener Ratio [-]");
ylim([1-1e-4,1+1e-4])
xlabel("$\phi$ [-]",'Interpreter','latex');

fontsize(gcf,40,'points')
lgd = legend('Hexagon','Reinforced hexagon');
lgd.FontSize = 30;

%% Figure 7 Homogenized bulk and shear
bulkHexa   = @(phi) (C11hexa(phi)-C33hexa(phi))./(C11hexa(0)-C33hexa(0));
shearHexa  = @(phi) C33hexa(phi)./(C33hexa(0));
bulkHoney  = @(phi) (C11honey(phi)-C33honey(phi))./(C11honey(0)-C33honey(0));
shearHoney = @(phi) C33honey(phi)./(C33honey(0));

bulkHexaMat = C11hexaMat - C33hexaMat;
shearHexaMat = C33hexaMat;
bulkHoneyMat = C11honeyMat - C33honeyMat;
shearHoneyMat = C33honeyMat;

figure(7)
t = tiledlayout(1,2);
nexttile
hold on
grid minor
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(bulkHexa,[0,1],'Color',cmp(4,:),'LineStyle','-','LineWidth',3);
fplot(bulkHoney,[0 1],'Color',cmp(4,:),'LineStyle','--','LineWidth',3);
ylabel('$\kappa(\phi)/\kappa_0$ [-]','Interpreter','latex');
ylim([0,inf])
xlabel("$\phi$ [-]",'Interpreter','latex');

nexttile
hold on
grid minor
p3 = fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p1 = fplot(shearHexa,[0 1],'Color',cmp(4,:),'LineStyle','-','LineWidth',3);
p2 = fplot(shearHoney,[0 1],'Color',cmp(4,:),'LineStyle','--','LineWidth',3);
ylabel('$\mu(\phi)/\mu_0$ [-]','Interpreter','latex');
ylim([0,inf])
xlabel("$\phi$ [-]",'Interpreter','latex');

fontsize(gcf,40,'points')
lgd = legend([p1,p2,p3],'Hexagon','Reinforced hexagon','H-S bound');
lgd.FontSize = 30;


%% Figure 7 initial derivative for different Poisson
kPrimeUB = @(nu) -1-k1(nu)/etak1(nu);
muPrimeUB = @(nu) -1-mu1(nu)/etamu1(nu);
gPrimeAT1 = @(nu) -2;
gPrimeAT2 = @(nu) -24;
gPrimeRational = @(nu,gamma) -gamma;

[dataHexa]  = load('HexagonDerivativeNu');
C11hexaNu = squeeze(dataHexa.mat(1,1,1,1,:));
C12hexaNu = squeeze(dataHexa.mat(2,2,1,1,:));
C33hexaNu = squeeze(dataHexa.mat(1,2,1,2,:));
[dataHoney] = load('HoneycombDerivativeNu');
C11honeyNu= squeeze(dataHoney.mat(1,1,1,1,:));
C12honeyNu= squeeze(dataHoney.mat(2,2,1,1,:));
C33honeyNu= squeeze(dataHoney.mat(1,2,1,2,:));
nuV = linspace(-0.99,0.5,200);

bulkHexaNu   = (C11hexaNu - C33hexaNu)./k1(nuV)';
shearHexaNu  = C33hexaNu./mu1(nuV)';
bulkHoneyNu  = (C11honeyNu - C33honeyNu)./k1(nuV)';
shearHoneyNu = C33honeyNu./mu1(nuV)';

bulkPrimeHexa   = -(1-bulkHexaNu)/dataHexa.phi;
shearPrimeHexa  = -(1-shearHexaNu)/dataHexa.phi;
bulkPrimeHoney  = -(1-bulkHoneyNu)/dataHoney.phi;
shearPrimeHoney = -(1-shearHoneyNu)/dataHoney.phi;

figure(8)
tiledlayout(1,2)
nexttile
hold on
grid minor
fplot(kPrimeUB,[-0.99,0.5],'Color',cmpGrad(1,:),'LineWidth',3)
fplot(gPrimeAT1,[-1,0.5],'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
fplot(gPrimeAT2,[-1,0.5],'Color',cmp(1,:),'LineStyle','--','LineWidth',3)
fplot(@(nu) gPrimeRational(nu,1),[-1,0.5],'Color',cmpGamma(6,:),'LineWidth',3)
fplot(@(nu) gPrimeRational(nu,3),[-1,0.5],'Color',cmpGamma(7,:),'LineWidth',3)
fplot(@(nu) gPrimeRational(nu,5),[-1,0.5],'Color',cmpGamma(8,:),'LineWidth',3)
fplot(@(nu) gPrimeRational(nu,9),[-1,0.5],'Color',cmpGamma(9,:),'LineWidth',3)
plot(nuV,bulkPrimeHexa,'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
plot(nuV,bulkPrimeHoney,'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel('$\kappa''(\phi)/\kappa_0$ [-]','Interpreter','latex');
xlabel("$\nu$ [-]",'Interpreter','latex');
ylim([-25 0])

nexttile
hold on
grid minor
p6 = fplot(muPrimeUB,[-0.81,0.5],'Color',cmpGrad(1,:),'LineWidth',3);
p1 = fplot(gPrimeAT1,[-1,0.5],'Color',cmp(1,:),'LineStyle','-','LineWidth',3);
p2 = fplot(gPrimeAT2,[-1,0.5],'Color',cmp(1,:),'LineStyle','--','LineWidth',3);
p3 = fplot(@(nu) gPrimeRational(nu,1),[-1,0.5],'Color',cmpGamma(6,:),'LineWidth',3);
fplot(@(nu) gPrimeRational(nu,3),[-1,0.5],'Color',cmpGamma(7,:),'LineWidth',3);
fplot(@(nu) gPrimeRational(nu,5),[-1,0.5],'Color',cmpGamma(8,:),'LineWidth',3);
fplot(@(nu) gPrimeRational(nu,9),[-1,0.5],'Color',cmpGamma(9,:),'LineWidth',3);
p4 = plot(nuV(24:end),shearPrimeHexa(24:end),'Color',cmp(4,:),'LineStyle','-','LineWidth',3);
p5 = plot(nuV(21:end),shearPrimeHoney(21:end),'Color',cmp(4,:),'LineStyle','--','LineWidth',3);
ylabel('$\mu''(\phi)/\mu_0$ [-]','Interpreter','latex');
xlabel("$\nu$ [-]",'Interpreter','latex');
ylim([-25 0])

fontsize(gcf,40,'points')
lgd = legend([p1,p2,p3,p4,p5,p6],'AT1','AT2','Rational','Hexagon','Reinforced hexagon','H-S Upper bound');
lgd.FontSize = 30;

eq = @(nu) kPrimeUB(nu)-muPrimeUB(nu);
fzero(eq,0.3)
kPrimeUB(fzero(eq,0.3))
muPrimeUB(fzero(eq,0.3))

%% Figure 9 Results 1Elem 
% (LOAD DATA)
resAT1 = load('1ElemAT1_v2.mat');
resAT2 = load('1ElemAT2_v2.mat');
resRat1 = load('1ElemRational198_v2.mat');
resRat2 = load('1ElemRational1_v2');
resHexa = load('1ElemHexagon_v2.mat');
resHoney = load('1ElemHoneycomb_v2.mat');

% (PLOT)

figure(9)
t = tiledlayout(2,1);
nexttile
hold on
grid minor
plot(resAT1.outputData.displacement.value,resAT1.outputData.force,'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
plot(resAT2.outputData.displacement.value,resAT2.outputData.force,'Color',cmp(1,:),'LineStyle','--','LineWidth',3)
plot([resRat1.outputData.displacement.value, 0.04],[resRat1.outputData.force,0],'Color',cmpGamma(6,:),'LineStyle','-','LineWidth',3)
plot([resRat2.outputData.displacement.value, 0.04],[resRat2.outputData.force,0],'Color',cmpGamma(6,:),'LineStyle','--','LineWidth',3)
plot(resHexa.outputData.displacement.value,resHexa.outputData.force,'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
plot(resHoney.outputData.displacement.value,resHoney.outputData.force,'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel("Force [kN]",'Interpreter','latex');
xlabel({"Displacement [mm]";"(a)"},'interpreter','latex')
xlim([0,0.04])

nexttile
hold on
grid minor
plot(resAT1.outputData.displacement.value,resAT1.outputData.damage.maxValue,'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
plot(resAT2.outputData.displacement.value,resAT2.outputData.damage.maxValue,'Color',cmp(1,:),'LineStyle','--','LineWidth',3)
plot([resRat1.outputData.displacement.value,0.04],[resRat1.outputData.damage.maxValue,1],'Color',cmpGamma(6,:),'LineStyle','-','LineWidth',3)
plot([resRat2.outputData.displacement.value,0.04],[resRat2.outputData.damage.maxValue,1],'Color',cmpGamma(6,:),'LineStyle','--','LineWidth',3)
plot(resHexa.outputData.displacement.value,resHexa.outputData.damage.maxValue,'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
plot(resHoney.outputData.displacement.value,resHoney.outputData.damage.maxValue,'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel("$\phi$ [-]",'Interpreter','latex');
xlabel({"Displacement [mm]";"(b)"},'interpreter','latex')
xlim([0,0.04])

fontsize(gcf,40,'points')
lgd = legend({'AT1','AT2','Rational ($\sigma_c$ = 1.984 GPa)','Rational ($\sigma_c$ = 1 GPa)','Hexagon','Reinforced Hexagon'},'Interpreter','latex','Location','eastoutside');
lgd.Layout.Tile = 'east';
lgd.FontSize = 30;

%% Figure 10 Results others

% resAT1 = load('UniaxialNewMeshAT1.mat');
% resAT2 = load('UniaxialNewMeshAT2.mat');
% resRat1 = load('UniaxialNewMeshRational3.mat');
% resRat2 = load('UniaxialNewMeshRational2.mat');
% resHexa = load('UniaxialNewMeshHexagon.mat');
% resHoney = load('UniaxialNewMeshHoneycomb.mat');

% resAT1 = load('SENtractionNewMeshAT1.mat');
% resAT2 = load('SENtractionNewMeshAT2.mat');
% resRat1 = load('SENtractionNewMeshRational24.mat');
% resRat2 = load('SENtractionNewMeshRational5.mat');
% resHexa = load('SENtractionNewMeshHexagon.mat');
% resHoney = load('SENtractionNewMeshHoneycomb.mat');

resAT1 = load('SENshearNewMeshAT1.mat');
resAT2 = load('SENshearNewMeshAT2.mat');
resRat1 = load('SENshearNewMeshRational24.mat');
resRat2 = load('SENshearNewMeshRational5.mat');
resHexa = load('SENshearNewMeshHexagon.mat');
resHoney = load('SENshearNewMeshHoneycomb.mat');


figure(10)
t = tiledlayout(1,2);
nexttile
hold on
grid minor
plot(resAT1.outputData.displacement.value,resAT1.outputData.force,'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
plot(resAT2.outputData.displacement.value,resAT2.outputData.force,'Color',cmp(1,:),'LineStyle','--','LineWidth',3)
plot([resRat1.outputData.displacement.value, 0.04],[resRat1.outputData.force,0],'Color',cmpGamma(6,:),'LineStyle','-','LineWidth',3)
plot([resRat2.outputData.displacement.value, 0.04],[resRat2.outputData.force,0],'Color',cmpGamma(6,:),'LineStyle','--','LineWidth',3)
plot(resHexa.outputData.displacement.value,resHexa.outputData.force,'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
plot(resHoney.outputData.displacement.value,resHoney.outputData.force,'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel("Force [kN]",'Interpreter','latex');
xlabel({"Displacement [mm]";"(a)"},'interpreter','latex')
%xlim([0,0.03])
%xlim([0,0.007])
xlim([0,0.02])
ax = gca;
ax.XAxis.Exponent = 0;

nexttile
hold on
grid minor
plot(resAT1.outputData.displacement.value,resAT1.outputData.damage.maxValue,'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
plot(resAT2.outputData.displacement.value,resAT2.outputData.damage.maxValue,'Color',cmp(1,:),'LineStyle','--','LineWidth',3)
plot([resRat1.outputData.displacement.value,0.04],[resRat1.outputData.damage.maxValue,1],'Color',cmpGamma(6,:),'LineStyle','-','LineWidth',3)
plot([resRat2.outputData.displacement.value,0.04],[resRat2.outputData.damage.maxValue,1],'Color',cmpGamma(6,:),'LineStyle','--','LineWidth',3)
plot(resHexa.outputData.displacement.value,resHexa.outputData.damage.maxValue,'Color',cmp(4,:),'LineStyle','-','LineWidth',3)
plot(resHoney.outputData.displacement.value,resHoney.outputData.damage.maxValue,'Color',cmp(4,:),'LineStyle','--','LineWidth',3)
ylabel("$\phi$ [-]",'Interpreter','latex');
xlabel({"Displacement [mm]";"(b)"},'interpreter','latex')
%xlim([0,0.03])
%xlim([0,0.007])
xlim([0,0.02])
ax = gca;
ax.XAxis.Exponent = 0;

fontsize(gcf,40,'points')
lgd = legend({'AT1','AT2','Rational ($\sigma_c$ = 2.445 GPa)','Rational ($\sigma_c$ = 5 GPa)','Hexagon','Reinforced Hexagon'},'Interpreter','latex','Location','eastoutside');
lgd.Layout.Tile = 'east';
lgd.FontSize = 30;



%% Prints

dmgAT1 = resAT1.outputData.damage.field.fun;
dmgAT2 = resAT2.outputData.damage.field.fun;
dmgRat1 = resRat1.outputData.damage.field.fun;
dmgRat2 = resRat2.outputData.damage.field.fun;
dmgHexa = resHexa.outputData.damage.field.fun;
dmgHoney = resHoney.outputData.damage.field.fun;


uAT1 = resAT1.outputData.displacement.field;
uAT2 = resAT2.outputData.displacement.field;
uRat1 = resRat1.outputData.displacement.field;
uRat2 = resRat2.outputData.displacement.field;
uHexa = resHexa.outputData.displacement.field;
uHoney = resHoney.outputData.displacement.field;

printResult(dmgAT1,uAT1,'AT1')
printResult(dmgAT2,uAT2,'AT2')
printResult(dmgRat1,uRat1,'Rat1')
printResult(dmgRat2,uRat2,'Rat2')
printResult(dmgHexa,uHexa,'Hexa')
printResult(dmgHoney,uHoney,'Honey')

function printResult(dmgFun,uFun,filename)
    s.mesh = dmgFun.mesh;
    s.fun = {dmgFun};
    s.type = 'Paraview';
    s.filename = filename;
    p = FunctionPrinter.create(s);
    p.appendFunction(uFun,'disp')
    p.print()
end