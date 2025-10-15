%% DEGRADATION FUNCTION
E = 210;
load('CircleAreaDerivative10.mat')
C10 = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('CircleAreaDerivative15.mat')
C15 = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('CircleAreaDerivative20.mat')
C20 = @(phi) E.*degradation.fun{1,1,1,1}(phi);

cmp = orderedcolors("gem");
figure()
hold on
fplot(C10,[0 1],'Color',cmp(1,:))
fplot(C15,[0 1],'Color',cmp(2,:))
fplot(C20,[0 1],'Color',cmp(3,:));
legend('$\sigma = 1MPa$','$\sigma = 1.5MPa$','$\sigma = 2MPa$','Interpreter','latex')
title('Circle Degradation function ($\nu = 0.3$)','Interpreter','latex')

%% 1 ELEMENT COMPARISON (sigma max)

load('Cir10.mat')
u10 = outputData.displacement.value;
f10 = outputData.force;
d10 = outputData.damage.maxValue;
load('Cir15.mat')
u15 = outputData.displacement.value;
f15 = outputData.force;
d15 = outputData.damage.maxValue;
load('Cir20.mat')
u20 = outputData.displacement.value;
f20 = outputData.force;
d20 = outputData.damage.maxValue;

cmp = orderedcolors("gem");
figure()
hold on
plot(u10,f10,'Color',cmp(1,:))
plot(u15,f15,'Color',cmp(2,:))
plot(u20,f20,'Color',cmp(3,:));
legend('$\sigma = 1MPa$','$\sigma = 1.5MPa$','$\sigma = 2MPa$','Interpreter','latex')
title('Circle One Element (Force-displacement)','Interpreter','latex')

figure()
hold on
plot(d10,f10,'Color',cmp(1,:))
plot(d15,f15,'Color',cmp(2,:))
plot(d20,f20,'Color',cmp(3,:));
legend('$\sigma = 1MPa$','$\sigma = 1.5MPa$','$\sigma = 2MPa$','Interpreter','latex')
title('Circle One Element (Force-damage)','Interpreter','latex')

figure()
hold on
plot(u10,d10,'Color',cmp(1,:))
plot(u15,d15,'Color',cmp(2,:))
plot(u20,d20,'Color',cmp(3,:));
legend('$\sigma = 1MPa$','$\sigma = 1.5MPa$','$\sigma = 2MPa$','Interpreter','latex')
title('Circle One Element (Damage)','Interpreter','latex')

%% 1 ELEMENT COMPARISON (poisson)

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