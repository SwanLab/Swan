%% Define SIMP-ALL and H-S Bounds formulas
close all

E  = 210;
nu = 0.3;
k  = E./(2.*(1-nu));
mu = E./(2.*(1+nu));

k0     = 1e-10;
k1     = k;
mu0    = 1e-10;
mu1    = mu;
etak0  = mu0;
etak1  = mu1;
etamu0 = (k0.*mu0)./(2.*mu0+k0);
etamu1 = (k1.*mu1)./(2.*mu1+k1);

muSIMPALL = @(phi) ((mu1-mu0).*(etamu0-etamu1).*(phi).*(1-phi) + mu0.*(mu1+etamu0).*(phi) + mu1.*(mu0+etamu1).*(1-phi))./...
                      ((mu1+etamu0).*(phi) + (mu0+etamu1).*(1-phi));

kSIMPALL = @(phi) ((k1-k0).*(etak0-etak1).*(phi).*(1-phi) + k0.*(k1+etak0).*(phi) + k1.*(k0+etak1).*(1-phi))./...
                         ((k1+etak0).*(phi) + (k0+etak1).*(1-phi));

kUB  = @(phi) k0.*(phi) + k1.*(1-phi) - ((1-phi).*phi.*(k1-k0).^2)./(k0.*(1-phi) + k1.*phi + etak1);
muUB = @(phi) mu0.*(phi) + mu1.*(1-phi) - ((1-phi).*phi.*(mu1-mu0).^2)./(mu0.*(1-phi) + mu1.*phi + etamu1);
kLB  = @(phi) k0.*(phi) + k1.*(1-phi) - ((1-phi).*phi.*(k1-k0).^2)./(k0.*(1-phi) + k1.*phi + etak0);
muLB = @(phi) mu0.*(phi) + mu1.*(1-phi) - ((1-phi).*phi.*(mu1-mu0).^2)./(mu0.*(1-phi) + mu1.*phi + etamu0);

%% Degradation functions
% AT degradation
kAT1  = @(phi) k.*(1-phi).^2;
muAT1 = @(phi) mu.*(1-phi).^2;

kAT2  = @(phi) k.*(1-sqrt(phi)).^2;
muAT2 = @(phi) mu.*(1-sqrt(phi)).^2;

% Rational (Wu)
a1 = -6;
kRatio  = @(phi) k.*(phi.^2 - 2.*phi + 1)./(1-(a1+2).*phi);
muRatio = @(phi) mu.*(phi.^2 - 2.*phi + 1)./(1-(a1+2).*phi);

% Rational (Alessi)
gammaK  = 5;
gammaMu = 1.5;
w = @(phi) 1 - (1-phi).^2;
kGamma = @(phi)  k.*(1-w(phi))/(1+(gammaK-1)*w(phi));
muGamma = @(phi) mu.*(1-w(phi))/(1+(gammaMu-1)*w(phi));

% Homogenized
load("HexagonArea.mat")

% C11 = squeeze(mat(1,1,1,1,:));
% C12 = squeeze(mat(1,1,2,2,:));
% C33 = squeeze(mat(1,2,1,2,:));
% ar = (C33)./(C11-C12);
% plot(ar)
% ylim([0.99999, 1.00001])

muHomog = @(phi) (210/2).*degradation.fun{1,2,1,2,:}(phi);
kHomog  = @(phi) (210).*degradation.fun{1,1,1,1,:}(phi) - muHomog(phi);


%% Plot
cmp = orderedcolors("gem");

figure(1)
hold on
fplot(kAT1,[0 1]);
fplot(kAT2,[0 1]);
fplot(kRatio,[0 1]);
fplot(kGamma,[0 1]);
fplot(kSIMPALL,[0 1]);
fplot(kHomog,[0 1])
legend('AT1','AT2','Wu','Alessi','SIMPALL','Hexagon')
title('Bulk modulus degradation')
ylabel('$\kappa$ [GPa]','Interpreter','latex')
xlabel('Damage $\phi$ [-]','Interpreter','latex')

figure(2)
hold on
fplot(kAT1,[0 1]);
fplot(kAT2,[0 1]);
fplot(kRatio,[0 1]);
fplot(kGamma,[0 1]);
fplot(kSIMPALL,[0 1]);
fplot(kHomog,[0 1])
legend('AT1','AT2','Wu','Alessi','SIMPALL','Hexagon')
title('Shear modulus degradation')
ylabel('$\mu$ [GPa]','Interpreter','latex')
xlabel('Damage $\phi$ [-]','Interpreter','latex')


%% Plot derivatives
syms phi

figure(3)
hold on
fplot(diff(kAT1(phi)),[0 1]);
fplot(diff(kAT2(phi)),[0 1]);
fplot(diff(kRatio(phi)),[0 1]);
fplot(diff(kGamma(phi)),[0 1]);
fplot(diff(kSIMPALL(phi)),[0 1]);
fplot(diff(kHomog(phi)),[0 1])
legend('AT1','AT2','Wu','Alessi','SIMPALL','Hexagon')
title('Bulk modulus degradation derivative')
ylabel('$\partial\mkappa / \partial\phi$ [GPa]','Interpreter','latex')
xlabel('Damage $\phi$ [-]','Interpreter','latex')

figure(4)
hold on
fplot(diff(muAT1(phi)),[0 1]);
fplot(diff(muAT2(phi)),[0 1]);
fplot(diff(muRatio(phi)),[0 1]);
fplot(diff(muGamma(phi)),[0 1]);
fplot(diff(muSIMPALL(phi)),[0 1]);
fplot(diff(muHomog(phi)),[0 1])
legend('AT1','AT2','Wu','Alessi','SIMPALL','Hexagon')
title('Shear modulus degradation derivative')
ylabel('$\partial\mu / \partial\phi$ [GPa]','Interpreter','latex')
xlabel('Damage $\phi$ [-]','Interpreter','latex')

figure(5)
hold on
fplot(diff(diff(kAT1(phi))),[0 1]);
fplot(diff(diff(kAT2(phi))),[0 1]);
fplot(diff(diff(kRatio(phi))),[0 1]);
fplot(diff(diff(kGamma(phi))),[0 1]);
fplot(diff(diff(kSIMPALL(phi))),[0 1]);
fplot(diff(diff(kHomog(phi))),[0 1])
legend('AT1','AT2','Wu','Alessi','SIMPALL','Hexagon')
title('Bulk modulus degradation second derivative')
ylabel('$\partial^2\mkappa / \partial\phi^2$ [GPa]','Interpreter','latex')
xlabel('Damage $\phi$ [-]','Interpreter','latex')

figure(6)
hold on
fplot(diff(diff(muAT1(phi))),[0 1]);
fplot(diff(diff(muAT2(phi))),[0 1]);
fplot(diff(diff(muRatio(phi))),[0 1]);
fplot(diff(diff(muGamma(phi))),[0 1]);
fplot(diff(diff(muSIMPALL(phi))),[0 1]);
fplot(diff(diff(muHomog(phi))),[0 1])
legend('AT1','AT2','Wu','Alessi','SIMPALL','Hexagon')
title('Shear modulus degradation second derivative')
ylabel('$\partial^2\mu / \partial\phi^2$ [GPa]','Interpreter','latex')
xlabel('Damage $\phi$ [-]','Interpreter','latex')