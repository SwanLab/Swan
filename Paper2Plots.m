close all
E  = 210;
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
figure(4)
tiledlayout(1,2)
nexttile
hold on
fplot(@(phi) gRational(phi,1),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) gRational(phi,3),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gRational(phi,5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gRational(phi,10),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gRational(phi,50),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
fontsize(gcf,25,'points')
title('Bulk modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$ [-]','Interpreter','latex');
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
nexttile
hold on
fplot(@(phi) gRational(phi,1),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) muUB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) gRational(phi,3),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gRational(phi,5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gRational(phi,10),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gRational(phi,50),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) muLB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
fontsize(gcf,25,'points')
title('Shear modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$ [-]','Interpreter','latex');
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
legend({'Wu','H-S bounds'},'Interpreter','latex')

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
tiledlayout(1,3)
nexttile
hold on
plot(dataHexa.phi,C11hexa)
plot(dataHoney.phi,C11honey)
ylabel(char(8450)+"11 [GPa]");
ylim([0,inf])
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
fontsize(gcf,25,'points')
nexttile
hold on
plot(dataHexa.phi,C12hexa)
plot(dataHoney.phi,C12honey)
ylabel(char(8450)+"12 [GPa]");
ylim([0,inf])
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
fontsize(gcf,25,'points')
nexttile
hold on
plot(dataHexa.phi,C33hexa)
plot(dataHoney.phi,C33honey)
ylabel(char(8450)+"33 [GPa]");
ylim([0,inf])
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
fontsize(gcf,25,'points')
legend('Hexagon','Reinforced hexagon')

%% Figure 6 Zener Ratio
ZenerRatioHexa = 2*C33hexa./(C11hexa-C12hexa);
ZenerRatioHoney = 2*C33honey./(C11honey-C12honey);

figure(6)
hold on
plot(dataHexa.phi,ZenerRatioHexa);
plot(dataHoney.phi,ZenerRatioHoney);
ylabel("Zener Ratio [-]");
%ylim([1-1e-5,1+1e-5])
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
fontsize(gcf,25,'points')
legend('Hexagon','Reinforced hexagon')

%% Figure 7 Homogenized bulk and shear
bulkHexa   = C11hexa - C33hexa;
shearHexa  = C33hexa;
bulkHoney  = C11honey - C33honey;
shearHoney = C33honey;

figure(7)
tiledlayout(1,2)
nexttile
hold on
plot(dataHexa.phi,bulkHexa)
plot(dataHoney.phi,bulkHoney)
ylabel("Bulk modulus ($\kappa$) [MPa]",'Interpreter','latex');
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
fontsize(gcf,25,'points')
nexttile
hold on
plot(dataHexa.phi,shearHexa)
plot(dataHoney.phi,shearHoney)
ylabel("Shear modulus ($\kappa$) [MPa]",'Interpreter','latex');
xlabel("Damage ($\phi$) [-]",'Interpreter','latex');
fontsize(gcf,25,'points')
legend('Hexagon','Reinforced hexagon')