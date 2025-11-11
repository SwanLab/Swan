close all
E  = 210;
nu = -0.5;
k  = E./(2.*(1-nu));
mu = E./(2.*(1+nu));
C11 = E/((1+nu)*(1-nu));

k0     = 1e-10;
k1     = k;
mu0    = 1e-10;
mu1    = mu;
etak0  = mu0;
etak1  = mu1;
etamu0 = (k0.*mu0)./(2.*mu0+k0);
etamu1 = (k1.*mu1)./(2.*mu1+k1);

Gc = 5e-3;
sigma = 1.5;
lch   = (2*E*Gc)/(sigma^2);
slope = (- k1 - k1^2/etak1 - mu1 - mu1^2/etamu1)/C11;
lhs   = -2*(3/8)*(Gc/slope)*E*(1/sigma)^2;
lmax = (3/8)*lch;
derivFactorLmax = -2*(3/8)*(Gc/lmax)*E*(1/sigma)^2;
derivFactorLHS = -2*(3/8)*(Gc/lhs)*E*(1/sigma)^2;
derivFactor05LHS = -2*(3/8)*(Gc/(0.5*lhs))*E*(1/sigma)^2;
derivFactor001LHS = -2*(3/8)*(Gc/(0.01*lhs))*E*(1/sigma)^2;


degSIMPmu = @(phi) ((mu1-mu0).*(etamu0-etamu1).*(phi).*(1-phi) + mu0.*(mu1+etamu0).*(phi) + mu1.*(mu0+etamu1).*(1-phi))./...
                      ((mu1+etamu0).*(phi) + (mu0+etamu1).*(1-phi));

degSIMPkappa = @(phi) ((k1-k0).*(etak0-etak1).*(phi).*(1-phi) + k0.*(k1+etak0).*(phi) + k1.*(k0+etak1).*(1-phi))./...
                         ((k1+etak0).*(phi) + (k0+etak1).*(1-phi));

kUB  = @(phi) k0.*(phi) + k1.*(1-phi) - ((1-phi).*phi.*(k1-k0).^2)./(k0.*(1-phi) + k1.*phi + etak1);
muUB = @(phi) mu0.*(phi) + mu1.*(1-phi) - ((1-phi).*phi.*(mu1-mu0).^2)./(mu0.*(1-phi) + mu1.*phi + etamu1);
kLB  = @(phi) k0.*(phi) + k1.*(1-phi) - ((1-phi).*phi.*(k1-k0).^2)./(k0.*(1-phi) + k1.*phi + etak0);
muLB = @(phi) mu0.*(phi) + mu1.*(1-phi) - ((1-phi).*phi.*(mu1-mu0).^2)./(mu0.*(1-phi) + mu1.*phi + etamu0);

%% Degradation functions

C11SIMP = @(phi) degSIMPkappa(phi) + degSIMPmu(phi);
C11UB   = @(phi) kUB(phi) + muUB(phi);
C11LB   = @(phi) kLB(phi) + muLB(phi);

C11AT1 = @(phi) C11.*(1-phi).^2; 
kAT1  = @(phi) k.*(1-phi).^2;
muAT1 = @(phi) mu.*(1-phi).^2;

C11AT2 = @(phi) C11.*(1-sqrt(phi)).^2;
kAT2  = @(phi) k.*(1-sqrt(phi)).^2;
muAT2 = @(phi) mu.*(1-sqrt(phi)).^2;

C11Linear = @(phi) C11.*(1-phi);
kLinear  = @(phi) k.*(1-phi);
muLinear = @(phi) mu.*(1-phi);

C11RationalLmax = @(phi) C11.*(phi.^2 - 2.*phi + 1)./(1-(derivFactorLmax+2).*phi);
kRational1  = @(phi) k.*(phi.^2 - 2.*phi + 1)./(1-(derivFactorLmax+2).*phi);
muRational1 = @(phi) mu.*(phi.^2 - 2.*phi + 1)./(1-(derivFactorLmax+2).*phi);

C11RationalLHS = @(phi) C11.*(phi.^2 - 2.*phi + 1)./(1-(derivFactorLHS+2).*phi);
kRational2  = @(phi) k.*(phi.^2 - 2.*phi + 1)./(1-(derivFactorLHS+2).*phi);
muRational2 = @(phi) mu.*(phi.^2 - 2.*phi + 1)./(1-(derivFactorLHS+2).*phi);

C11Rational05LHS = @(phi) C11.*(phi.^2 - 2.*phi + 1)./(1-(derivFactor05LHS+2).*phi);
C11Rational001LHS = @(phi) C11.*(phi.^2 - 2.*phi + 1)./(1-(derivFactor001LHS+2).*phi);

load('DegSqr15lHS.mat')
C11square = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('CircleArea.mat')
C11circle = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('HexagonArea.mat')
C11hexa = @(phi) E.*degradation.fun{1,1,1,1}(phi);

%% Plot degradation C11
cmp = orderedcolors("gem");

figure()
hold on
fplot(C11AT1,[0 1],'Color',cmp(1,:))
fplot(C11AT2,[0 1],'Color',cmp(1,:),'LineStyle','--')
fplot(C11Linear,[0 1],'Color',cmp(2,:));
fplot(C11RationalLmax,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker',"o")
fplot(C11RationalLHS,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker','*')
fplot(C11Rational05LHS,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker','+')
fplot(C11Rational001LHS,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker','x')
fplot(C11SIMP,[0 1],'Color',cmp(4,:))
fplot(C11UB,[0 1],'Color',cmp(4,:),'LineStyle','--')
fplot(C11LB,[0 1],'Color',cmp(4,:),'LineStyle','--')
fplot(C11square,[0 1])
legend('AT1','AT2','Linear','Rational (l0 = cw*lch = 0.35)','Rational (l0 = lHS)','Rational (l0 = 0.5lHS)','Rational (l0 = 0.01lHS)','Rational (2MPa)','SIMP','HS')
title('Degradation function ($\nu = 0.3$)','Interpreter','latex')

%% Plot degradation mu

figure()
hold on
fplot(muAT1,[0 1],'Color',cmp(1,:))
fplot(muAT2,[0 1],'Color',cmp(1,:),'LineStyle','--')
fplot(muLinear,[0 1],'Color',cmp(2,:));
fplot(muRational1,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker',"o")
fplot(muRational2,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker',"square")
fplot(degSIMPmu,[0 1],'Color',cmp(4,:))
fplot(muUB,[0 1],'Color',cmp(4,:),'LineStyle','--')
fplot(muLB,[0 1],'Color',cmp(4,:),'LineStyle','--')
legend('AT1','AT2','Linear','Rational (1MPa)','Rational (2MPa)','SIMP','HS')
title('Shear degradation function ($\nu = -0.5$)','Interpreter','latex')


%% Plot degradation kappa

figure()
hold on
fplot(kAT1,[0 1],'Color',cmp(1,:))
fplot(kAT2,[0 1],'Color',cmp(1,:),'LineStyle','--')
fplot(kLinear,[0 1],'Color',cmp(2,:));
fplot(kRational1,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker',"o")
fplot(kRational2,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker',"square")
fplot(degSIMPkappa,[0 1],'Color',cmp(4,:))
fplot(kUB,[0 1],'Color',cmp(4,:),'LineStyle','--')
fplot(kLB,[0 1],'Color',cmp(4,:),'LineStyle','--')
legend('AT1','AT2','Linear','Rational (1MPa)','Rational (2MPa)','SIMP','HS')
title('Bulk degradation function ($\nu = -0.5$)','Interpreter','latex')





