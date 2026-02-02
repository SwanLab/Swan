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

kUB  = @(phi,nu) (1/k1(nu))*(k0.*(phi) + k1(nu).*(1-phi) - ((1-phi).*phi.*(k1(nu)-k0).^2)./(k0.*(1-phi) + k1(nu).*phi + etak1(nu)));
muUB = @(phi,nu) (1/mu1(nu))*(mu0.*(phi) + mu1(nu).*(1-phi) - ((1-phi).*phi.*(mu1(nu)-mu0).^2)./(mu0.*(1-phi) + mu1(nu).*phi + etamu1(nu)));
kLB  = @(phi,nu) (1/k1(nu))*(k0.*(phi) + k1(nu).*(1-phi) - ((1-phi).*phi.*(k1(nu)-k0).^2)./(k0.*(1-phi) + k1(nu).*phi + etak0));
muLB = @(phi,nu) (1/mu1(nu))*(mu0.*(phi) + mu1(nu).*(1-phi) - ((1-phi).*phi.*(mu1(nu)-mu0).^2)./(mu0.*(1-phi) + mu1(nu).*phi + etamu0));

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
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','diamond','MarkerSize',7,'LineWidth',3)
fplot(gAT2,[0 1],'Color',cmp(2,:),'LineStyle','-','Marker','o','MarkerSize',7,'LineWidth',3)
fontsize(gcf,40,'points')
ylabel('$\kappa(\phi)/\kappa_0$','Interpreter','latex');
nexttile
hold on
grid minor

fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p3 = fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p1 = fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','diamond','MarkerSize',7,'LineWidth',3)
p2 = fplot(gAT2,[0 1],'Color',cmp(2,:),'LineStyle','-','Marker','o','MarkerSize',7,'LineWidth',3)
fontsize(gcf,40,'points')
ylabel('$\mu(\phi)/\mu_0$','Interpreter','latex');
xlabel(t,"Damage ($\phi$)",'Interpreter','latex');
legend([p1,p2,p3],'AT1','AT2','H-S bounds')

%% Figure 2: AT vs H-S bounds for different Possion

cmpGrad = gray(10)
figure(2)
t = tiledlayout(1,2);
nexttile
hold on
grid minor
fplot(@(phi) kUB(phi,-0.9),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,-0.5),[0 1],'Color',cmpGrad(2,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmpGrad(3,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0.3),[0 1],'Color',cmpGrad(4,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0.5),[0 1],'Color',cmpGrad(5,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(6,:),'LineStyle','-','LineWidth',3);
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','diamond','MarkerSize',7,'LineWidth',3)
fplot(gAT2,[0 1],'Color',cmp(2,:),'LineStyle','-','Marker','o','MarkerSize',7,'LineWidth',3)
fontsize(gcf,40,'points')
ylabel('$\kappa(\phi)/\kappa_0$','Interpreter','latex');
nexttile
hold on
grid minor
p3 = fplot(@(phi) muUB(phi,-0.9),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,-0.5),[0 1],'Color',cmpGrad(2,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,0),[0 1],'Color',cmpGrad(3,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmpGrad(4,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muUB(phi,0.5),[0 1],'Color',cmpGrad(5,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p1 = fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','diamond','MarkerSize',7,'LineWidth',3)
p2 = fplot(gAT2,[0 1],'Color',cmp(2,:),'LineStyle','-','Marker','o','MarkerSize',7,'LineWidth',3)
fontsize(gcf,40,'points')
ylabel('$\mu(\phi)/\mu_0$','Interpreter','latex');
xlabel(t,"Damage ($\phi$)",'Interpreter','latex');
legend([p1,p2,p3],'AT1','AT2','H-S bounds')

%% Alessi and Wu functions
gRational = @(phi,gamma) (1-phi)/((1-phi)+gamma*phi);


%% Figure 4_ Rational vs H-S bounds
cmpGamma = sky(10);

figure(4)
t = tiledlayout(1,2);
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
fontsize(gcf,40,'points')
ylabel('$\kappa(\phi)/\kappa_0$','Interpreter','latex');
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
fontsize(gcf,40,'points')
ylabel('$\mu(\phi)/\mu_0$','Interpreter','latex');
xlabel(t,"Damage ($\phi$)",'Interpreter','latex');
legend([p1,p2],{'Rational','H-S bounds'},'Interpreter','latex')

%% Figure 5 Constitutive tensor homogenized
[dataHexa]  = load('HexagonAreaNew.mat');
C11hexa = squeeze(dataHexa.mat(1,1,1,1,:));
C12hexa = squeeze(dataHexa.mat(2,2,1,1,:));
C33hexa = squeeze(dataHexa.mat(1,2,1,2,:));
[dataHoney] = load('HoneycombAreaNew.mat');
C11honey= squeeze(dataHoney.mat(1,1,1,1,:));
C12honey= squeeze(dataHoney.mat(2,2,1,1,:));
C33honey= squeeze(dataHoney.mat(1,2,1,2,:));

figure(5)
t = tiledlayout(1,3);
nexttile
hold on
grid minor
fplot(dataHexa.degradation.fun{1,1,1,1},[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
fplot(dataHoney.degradation.fun{1,1,1,1},[0 1],'Color',cmp(2,:),'LineStyle','-','LineWidth',3)
ylabel(char(8450)+"11 [GPa]");
ylim([0,inf])
fontsize(gcf,40,'points')
nexttile
hold on
grid minor
fplot(dataHexa.degradation.fun{2,2,1,1},[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
fplot(dataHoney.degradation.fun{2,2,1,1},[0 1],'Color',cmp(2,:),'LineStyle','-','LineWidth',3)
ylabel(char(8450)+"12 [GPa]");
ylim([0,inf])
fontsize(gcf,40,'points')
nexttile
hold on
grid minor
fplot(dataHexa.degradation.fun{1,2,1,2},[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
fplot(dataHoney.degradation.fun{1,2,1,2},[0 1],'Color',cmp(2,:),'LineStyle','-','LineWidth',3)
ylabel(char(8450)+"33 [GPa]");
ylim([0,inf])
xlabel(t,"Damage ($\phi$)",'Interpreter','latex');
fontsize(gcf,25,'points')
legend('Hexagon','Reinforced hexagon')

%% Figure 6 Zener Ratio
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
fplot(ZenerRatioHexa,[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3);
fplot(ZenerRatioHoney,[0 1],'Color',cmp(2,:),'LineStyle','-','LineWidth',3);
ylabel("Zener Ratio");
ylim([1-1e-4,1+1e-4])
xlabel("Damage ($\phi$)",'Interpreter','latex');
fontsize(gcf,40,'points')
legend('Hexagon','Reinforced hexagon')

%% Figure 7 Homogenized bulk and shear
bulkHexa   = @(phi) (C11hexa(phi)-C33hexa(phi))/(C11hexa(0)-C33hexa(0));
shearHexa  = @(phi) C33hexa(phi)/C33hexa(0);
bulkHoney  = @(phi) (C11honey(phi)-C33honey(phi))/(C11honey(0)-C33honey(0));
shearHoney = @(phi) C33honey(phi)/C33honey(0);

figure(7)
t = tiledlayout(1,2);
nexttile
hold on
grid minor
p3 = fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
p1 = fplot(bulkHexa,[0,1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3);
p2 = fplot(bulkHoney,[0 1],'Color',cmp(2,:),'LineStyle','-','LineWidth',3);
ylabel('$\kappa(\phi)/\kappa_0$','Interpreter','latex');
fontsize(gcf,40,'points')
nexttile
hold on
grid minor
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmpGrad(1,:),'LineStyle','-','LineWidth',3);
fplot(shearHexa,[0 1],'Color',cmp(1,:),'LineStyle','-','LineWidth',3)
fplot(shearHoney,[0 1],'Color',cmp(2,:),'LineStyle','-','LineWidth',3)
ylabel('$\mu(\phi)/\mu_0$','Interpreter','latex');
xlabel(t,"Damage ($\phi$)",'Interpreter','latex');
fontsize(gcf,40,'points')
legend([p1,p2,p3],'Hexagon','Reinforced hexagon','H-S bound')


%% Figure 8 initial derivative for different Poisson
kPrimeUB  = @(nu) (-k1(nu)-(k1(nu)^2)/etak1(nu))*(1/k1(nu));
muPrimeUB = @(nu) (-mu1(nu)-(mu1(nu)^2)/etamu1(nu))*(1/mu1(nu));
gPrimeAT1 = @(nu) -2;
gPrimeAT2 = @(nu) -Inf;
gPrimeRational = @(nu,gamma) -gamma;

[dataHexa]  = load('HexaNu.mat');
C11hexaNu = squeeze(dataHexa.mat(1,1,1,1,:));
C12hexaNu = squeeze(dataHexa.mat(2,2,1,1,:));
C33hexaNu = squeeze(dataHexa.mat(1,2,1,2,:));
[dataHoney] = load('HoneyNu.mat');
C11honeyNu= squeeze(dataHoney.mat(1,1,1,1,:));
C12honeyNu= squeeze(dataHoney.mat(2,2,1,1,:));
C33honeyNu= squeeze(dataHoney.mat(1,2,1,2,:));
nuV = linspace(-0.99,0.5,100);

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
fplot(kPrimeUB,[-0.99,0.5],'Color',cmp(4,:))
fplot(gPrimeAT1,[-1,0.5],'Color',cmp(1,:),'LineStyle','-','Marker','+')
fplot(gPrimeAT2,[-1,0.5],'Color',cmp(1,:),'LineStyle','-','Marker','o')
fplot(@(nu) gPrimeRational(nu,1),[-1,0.5],'Color',cmp(2,:))
plot(nuV,bulkPrimeHexa,'Color',cmp(3,:),'LineStyle','-','Marker','+')
plot(nuV,bulkPrimeHoney,'Color',cmp(3,:),'LineStyle','-','Marker','o')
fplot(@(nu) gPrimeRational(nu,3),[-1,0.5],'Color',cmp(2,:))
fplot(@(nu) gPrimeRational(nu,5),[-1,0.5],'Color',cmp(2,:))
fplot(@(nu) gPrimeRational(nu,10),[-1,0.5],'Color',cmp(2,:))
ylabel("Bulk modulus initial derivative $(g_{\kappa}'(0))$",'Interpreter','latex');
xlabel("Poisson ratio ($\nu$) [-]",'Interpreter','latex');
ylim([-25 0])
fontsize(gcf,40,'points')
nexttile
hold on
fplot(muPrimeUB,[-0.99,0.5],'Color',cmp(4,:))
fplot(gPrimeAT1,[-1,0.5],'Color',cmp(1,:),'LineStyle','-','Marker','+')
fplot(gPrimeAT2,[-1,0.5],'Color',cmp(1,:),'LineStyle','-','Marker','o')
fplot(@(nu) gPrimeRational(nu,1),[-1,0.5],'Color',cmp(2,:))
plot(nuV,shearPrimeHexa,'Color',cmp(3,:),'LineStyle','-','Marker','+')
plot(nuV,shearPrimeHoney,'Color',cmp(3,:),'LineStyle','-','Marker','o')
fplot(@(nu) gPrimeRational(nu,3),[-1,0.5],'Color',cmp(2,:))
fplot(@(nu) gPrimeRational(nu,5),[-1,0.5],'Color',cmp(2,:))
fplot(@(nu) gPrimeRational(nu,10),[-1,0.5],'Color',cmp(2,:))
ylim([-25 0])
ylabel("Shear modulus initial derivative ($(g_{\mu}'(0)$)",'Interpreter','latex');
xlabel("Poisson ratio ($\nu$) [-]",'Interpreter','latex');
fontsize(gcf,40,'points')
legend('H-S upper bound','AT1','AT2','Rational','Hexagon','Reinforced hexagon')

eq = @(nu) kPrimeUB(nu)-muPrimeUB(nu);
fzero(eq,0.3)
kPrimeUB(fzero(eq,0.3))
muPrimeUB(fzero(eq,0.3))