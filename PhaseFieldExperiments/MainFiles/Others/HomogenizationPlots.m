%% DEGRADATION FUNCTION
E = 210;
load('DegSqr15L0Lch.mat')
Lch = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('DegSqr15L0LHS.mat')
LHS = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('DegSqr15L02Lch.mat')
Lch2 = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('DegSqr15L02Lch.mat')
Lch095 = @(phi) E.*degradation.fun{1,1,1,1}(phi);

cmp = orderedcolors("gem");
figure()
hold on
fplot(Lch,[0 1],'Color',cmp(1,:))
fplot(Lch2,[0 1],'Color',cmp(2,:));
fplot(Lch095,[0 1],'Color',cmp(3,:));
fplot(LHS,[0 1],'Color',cmp(4,:))
legend('Lch','Lch2','Lch095','LHS','Interpreter','latex')
title('Square Degradation function ($\nu = 0.3$)','Interpreter','latex')

%% 1 ELEMENT COMPARISON (same type - sigma max)

load('RationalMaxL0.mat')
u10 = outputData.displacement.value;
f10 = outputData.force;
d10 = outputData.damage.maxValue;
load('RationalLHS.mat')
u15 = outputData.displacement.value;
f15 = outputData.force;
d15 = outputData.damage.maxValue;
load('Rational05LHS.mat')
u20 = outputData.displacement.value;
f20 = outputData.force;
d20 = outputData.damage.maxValue;
load('Rational001LHS.mat')
u30 = outputData.displacement.value;
f30 = outputData.force;
d30 = outputData.damage.maxValue;

cmp = orderedcolors("gem");
figure()
hold on
plot(u10,f10,'Color',cmp(1,:))
plot(u15,f15,'Color',cmp(2,:))
plot(u20,f20,'Color',cmp(3,:))
plot(u30,f30,'Color',cmp(4,:));
legend('l0 = cw*lch = 0.35','l0 = lHS','l0 = 0.5lHS','l0 = 0.01lHS')
title('Circle One Element (Force-displacement)','Interpreter','latex')

figure()
hold on
plot(d10,f10,'Color',cmp(1,:))
plot(d15,f15,'Color',cmp(2,:))
plot(d20,f20,'Color',cmp(3,:))
plot(d30,f30,'Color',cmp(4,:));
legend('l0 = cw*lch = 0.35','l0 = lHS','l0 = 0.5lHS','l0 = 0.01lHS')
title('Circle One Element (Force-damage)','Interpreter','latex')

figure()
hold on
plot(u10,d10,'Color',cmp(1,:))
plot(u15,d15,'Color',cmp(2,:))
plot(u20,d20,'Color',cmp(3,:))
plot(u30,d30,'Color',cmp(4,:));
legend('l0 = cw*lch = 0.35','l0 = lHS','l0 = 0.5lHS','l0 = 0.01lHS')
title('Circle One Element (Damage)','Interpreter','latex')


%% 1 ELEMENT COMPARISON (same type - poisson)

load('AT1nu05Clamped.mat')
u05 = outputData.displacement.value;
f05 = outputData.force;
d05 = outputData.damage.maxValue;
load('AT1nu03Clamped.mat')
u03 = outputData.displacement.value;
f03 = outputData.force;
d03 = outputData.damage.maxValue;
load('AT1nu0Clamped.mat')
u0 = outputData.displacement.value;
f0 = outputData.force;
d0 = outputData.damage.maxValue;
load('AT1nu_05Clamped.mat')
u_05 = outputData.displacement.value;
f_05 = outputData.force;
d_05 = outputData.damage.maxValue;

cmp = orderedcolors("gem");
figure()
hold on
plot(u05,f05,'Color',cmp(1,:))
plot(u03,f03,'Color',cmp(2,:))
plot(u0,f0,'Color',cmp(3,:))
plot(u_05,f_05,'Color',cmp(4,:))
legend('$\nu = 0.5$','$\nu = 0.3$','$\nu = 0$','$\nu = -0.5$','Interpreter','latex')
title('AT1 One Element (Force-displacement)','Interpreter','latex')

figure()
hold on
plot(d05,f05,'Color',cmp(1,:))
plot(d03,f03,'Color',cmp(2,:))
plot(d0,f0,'Color',cmp(3,:))
plot(d_05,f_05,'Color',cmp(4,:))
legend('$\nu = 0.5$','$\nu = 0.3$','$\nu = 0$','$\nu = -0.5$','Interpreter','latex')
title('AT1 One Element (Force-damage)','Interpreter','latex')

figure()
hold on
plot(u05,d05,'Color',cmp(1,:))
plot(u03,d03,'Color',cmp(2,:))
plot(u0,d0,'Color',cmp(3,:))
plot(u_05,d_05,'Color',cmp(4,:))
legend('$\nu = 0.5$','$\nu = 0.3$','$\nu = 0$','$\nu = -0.5$','Interpreter','latex')
title('AT1 One Element (damage-displacement)','Interpreter','latex')

%% 1 ELEMENT COMPARISON (all types - poisson)
load('AT1nu_05Clamped.mat')
uAT1 = outputData.displacement.value;
dAT1 = outputData.damage.maxValue;
fAT1 = outputData.force;
load('AT2nu_05Clamped.mat')
uAT2 = outputData.displacement.value;
dAT2 = outputData.damage.maxValue;
fAT2 = outputData.force;
load('Rat1nu_05Clamped.mat')
uRat1 = outputData.displacement.value;
dRat1 = outputData.damage.maxValue;
fRat1 = outputData.force;
load('Rat15nu_05Clamped.mat')
uRat15 = outputData.displacement.value;
dRat15 = outputData.damage.maxValue;
fRat15 = outputData.force;
load('Rat2nu_05Clamped.mat')
uRat2 = outputData.displacement.value;
dRat2 = outputData.damage.maxValue;
fRat2 = outputData.force;
load('SIMPALLnu_05Clamped.mat')
uSimp = outputData.displacement.value;
dSimp = outputData.damage.maxValue;
fSimp = outputData.force;

cmp = orderedcolors("gem");
figure()
hold on
plot(uAT1,fAT1,'Color',cmp(1,:))
plot(uAT2,fAT2,'Color',cmp(1,:),'LineStyle','--')
plot(uRat1,fRat1,'Color',cmp(3,:),"Marker","o","MarkerIndices",[1:2:length(fRat1)])
plot(uRat15,fRat15,'Color',cmp(3,:),"Marker","+","MarkerIndices",[1:2:length(fRat15)])
plot(uRat2,fRat2,'Color',cmp(3,:),"Marker","*","MarkerIndices",[1:2:length(fRat2)])
plot(uSimp,fSimp,'Color',cmp(4,:))
legend('AT1','AT2','Rational (1MPa)','Rational (1.5MPa)','Rational (2MPa)','SIMPALL')
title('Force-displacement 1Elem (nu = -0.5)')

figure()
hold on
plot(dAT1,fAT1,'Color',cmp(1,:))
plot(dAT2,fAT2,'Color',cmp(1,:),'LineStyle','--')
plot(dRat1,fRat1,'Color',cmp(3,:),"Marker","o","MarkerIndices",[1:2:length(fRat1)])
plot(dRat15,fRat15,'Color',cmp(3,:),"Marker","+","MarkerIndices",[1:2:length(fRat15)])
plot(dRat2,fRat2,'Color',cmp(3,:),"Marker","*","MarkerIndices",[1:2:length(fRat2)])
plot(dSimp,fSimp,'Color',cmp(4,:))
legend('AT1','AT2','Rational (1MPa)','Rational (1.5MPa)','Rational (2MPa)','SIMPALL')
title('Force-damage 1Elem (nu = -0.5)')

figure()
hold on
plot(uAT1,dAT1,'Color',cmp(1,:))
plot(uAT2,dAT2,'Color',cmp(1,:),'LineStyle','--')
plot(uRat1,dRat1,'Color',cmp(3,:),"Marker","o","MarkerIndices",[1:2:length(fRat1)])
plot(uRat15,dRat15,'Color',cmp(3,:),"Marker","+","MarkerIndices",[1:2:length(fRat15)])
plot(uRat2,dRat2,'Color',cmp(3,:),"Marker","*","MarkerIndices",[1:2:length(fRat2)])
plot(uSimp,dSimp,'Color',cmp(4,:))
legend('AT1','AT2','Rational (1MPa)','Rational (1.5MPa)','Rational (2MPa)','SIMPALL')
title('Damage-displacement 1Elem (nu = -0.5)')
