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

figure(1)
tiledlayout(1,2)
nexttile
hold on
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','+')
fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','o')
fplot(@(phi) kUB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Bulk modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
nexttile
hold on
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','+')
fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','o')
fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Shear modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
legend('AT1','AT2','H-S bounds')

%% Figure 2: AT vs H-S bounds for different Possion
figure(2)
tiledlayout(1,2)
nexttile
hold on
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','+')
fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','o')
fplot(@(phi) kUB(phi,-0.9),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) kUB(phi,-0.5),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) kUB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) kUB(phi,0.5),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) kUB(phi,0.9),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Bulk modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
nexttile
hold on
fplot(gAT1,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','+')
fplot(gAT2,[0 1],'Color',cmp(1,:),'LineStyle','-','Marker','o')
fplot(@(phi) muUB(phi,-0.9),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) muUB(phi,-0.5),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) muUB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) muUB(phi,0.5),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) muUB(phi,0.9),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Shear modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
legend('AT1','AT2','H-S bounds')

%% Alessi and Wu functions
gAlessi = @(phi,gamma) (1-phi)/(1+(gamma-1)*phi);
gWu = @(phi,gamma) (phi^2 - 2*phi + 1)/(1+(gamma-2)*phi);

%% Figure 3: Alessi vs H-S bounds 
figure(3)
tiledlayout(1,2)
nexttile
hold on
fplot(@(phi) gAlessi(phi,0.5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) kUB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) gAlessi(phi,1),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,3),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,10),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,50),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Bulk modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
nexttile
hold on
fplot(@(phi) gAlessi(phi,0.5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) muUB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) gAlessi(phi,1),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,3),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,10),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gAlessi(phi,50),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) muLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Shear modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
legend({'Alessi','H-S bounds'},'Interpreter','latex')

%% Figure 4_ Wu vs H-S bounds
figure(4)
tiledlayout(1,2)
nexttile
hold on
fplot(@(phi) gWu(phi,1),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) kUB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) gWu(phi,3),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gWu(phi,5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gWu(phi,10),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gWu(phi,50),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) kLB(phi,0.3),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Bulk modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
nexttile
hold on
fplot(@(phi) gWu(phi,1),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) muUB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
fplot(@(phi) gWu(phi,3),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gWu(phi,5),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gWu(phi,10),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) gWu(phi,50),[0 1],'Color',cmp(1,:),'LineStyle','-')
fplot(@(phi) muLB(phi,0),[0 1],'Color',cmp(4,:),'LineStyle','-');
title('Shear modulus degradation function ($\nu = 0.3$)','Interpreter','latex')
ylabel('Degradation $g(\phi)$','Interpreter','latex');
xlabel('Damage $\phi$','Interpreter','latex');
legend({'Wu','H-S bounds'},'Interpreter','latex')