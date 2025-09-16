E  = 210;
nu = 0.3;
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

Gc = 5e-3; l0 = 0.1;
derivFactor1 = -2*(3/8)*(Gc/l0)*E*(1/1)^2;
derivFactor2 = -2*(3/8)*(Gc/l0)*E*(1/2)^2;

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

C11AT = @(phi) C11.*(1-phi).^2;
C11Linear = @(phi) C11.*(1-phi);
C11Rational1 = @(phi) C11.*(phi.^2 - 2.*phi + 1)./(1-(derivFactor1+2).*phi);
C11Rational2 = @(phi) C11.*(phi.^2 - 2.*phi + 1)./(1-(derivFactor2+2).*phi);

load('SquarePerimeter.mat')
C11square = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('CirclePerimeter.mat')
C11circle = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('HexagonPerimeter.mat')
C11hexa = @(phi) E.*degradation.fun{1,1,1,1}(phi);

%% Plots
cmp = orderedcolors("gem");

figure()
hold on
fplot(C11AT,[0 1],'Color',cmp(1,:))
fplot(C11Linear,[0 1],'Color',cmp(2,:));
fplot(C11Rational1,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker',"o")
fplot(C11Rational2,[0 1],'Color',cmp(3,:),'LineStyle','--','Marker',"square")
fplot(C11SIMP,[0 1],'Color',cmp(4,:))
fplot(C11UB,[0 1],'Color',cmp(4,:),'LineStyle','--')
fplot(C11LB,[0 1],'Color',cmp(4,:),'LineStyle','--')
 fplot(C11square,[0 1],'Color',cmp(5,:),'LineStyle','--','Marker',"square")
 fplot(C11circle,[0 1],'Color',cmp(5,:),'LineStyle','--','Marker',"o")
 fplot(C11hexa,[0 1],'Color',cmp(5,:),'LineStyle','--','Marker',"hexagram")
%legend('AT','Linear','Rational (1MPa)','Rational (2MPa)','SIMP','HS')
legend('AT','Linear','Rational (1MPa)','Rational (2MPa)','SIMP','HS','Homogenized (Square)','Homogenized (Circle)','Homogenized (Hexagon)')
title('Degradation function ($\nu = 0.3$)','Interpreter','latex')


%% 1ELEM
figure()
plot(uAT,FAT,'Color',cmp(1,:))
hold on
plot(uLinear,Flinear,'Color',cmp(2,:))
plot(uRatio1,Fratio1,'Color',cmp(3,:),'LineStyle','--','Marker',"o","MarkerIndices",[1:2:length(Fratio1)])
plot(uRatio2,Fratio2,'Color',cmp(3,:),'LineStyle','--','Marker',"square","MarkerIndices",[1:5:length(Fratio2)])
plot(uSIMP,Fsimp,'Color',cmp(4,:))
plot(uSquare,Fsquare,'Color',cmp(5,:),'LineStyle','--','Marker',"square","MarkerIndices",[1:10:length(Fsquare)])
plot(uCircle,FCircle,'Color',cmp(5,:),'LineStyle','--','Marker',"o","MarkerIndices",[1:10:length(FCircle)])
plot(uHexa,Fhexa,'Color',cmp(5,:),'LineStyle','--','Marker',"hexagram","MarkerIndices",[1:10:length(Fhexa)])
legend('AT1','Linear','Rational (1MPa)','Rational (2MPa)','SIMP','Homogenized (Square)','Homogenized (Circle)','Homogenized (Hexagon)')
title('Force-displacement 1Elem')

figure()
plot(uAT,dAT,'Color',cmp(1,:))
hold on
plot(uLinear,dLinear,'Color',cmp(2,:))
plot(uRatio1,dRatio1,'Color',cmp(3,:),'LineStyle','--','Marker',"o","MarkerIndices",[1:2:length(Fratio1)])
plot(uRatio2,dRatio2,'Color',cmp(3,:),'LineStyle','--','Marker',"square","MarkerIndices",[1:5:length(Fratio2)])
plot(uSIMP,dSIMP,'Color',cmp(4,:))
plot(uSquare,dSquare,'Color',cmp(5,:),'LineStyle','--','Marker',"square","MarkerIndices",[1:10:length(Fsquare)])
plot(uCircle,dCircle,'Color',cmp(5,:),'LineStyle','--','Marker',"o","MarkerIndices",[1:10:length(FCircle)])
plot(uHexa,dHexa,'Color',cmp(5,:),'LineStyle','--','Marker',"hexagram","MarkerIndices",[1:10:length(Fhexa)])
legend('AT1','Linear','Rational (1MPa)','Rational (2MPa)','SIMP','Homogenized (Square)','Homogenized (Circle)','Homogenized (Hexagon)')
title('Damage-displacement 1Elem')





