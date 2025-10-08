E = 210;
load('SquareAreaDerivative10.mat')
C10 = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('SquareAreaDerivative15.mat')
C15 = @(phi) E.*degradation.fun{1,1,1,1}(phi);
load('SquareAreaDerivative20.mat')
C20 = @(phi) E.*degradation.fun{1,1,1,1}(phi);

cmp = orderedcolors("gem");
figure()
hold on
fplot(C10,[0 1],'Color',cmp(1,:))
fplot(C15,[0 1],'Color',cmp(2,:))
fplot(C20,[0 1],'Color',cmp(3,:));
legend('$\sigma = 1MPa$','$\sigma = 1.5MPa$','$\sigma = 2MPa$','Interpreter','latex')
title('Square Degradation function ($\nu = 0.3$)','Interpreter','latex')

load('Sqr10.mat')
u10 = outputData.displacement.value;
f10 = outputData.force;
d10 = outputData.damage.maxValue;
load('Sqr15.mat')
u15 = outputData.displacement.value;
f15 = outputData.force;
d15 = outputData.damage.maxValue;
load('Sqr20.mat')
u20 = outputData.displacement.value;
f20 = outputData.force;
d20 = outputData.damage.maxValue;

figure()
hold on
plot(u10,f10,'Color',cmp(1,:))
plot(u15,f15,'Color',cmp(2,:))
plot(u20,f20,'Color',cmp(3,:));
legend('$\sigma = 1MPa$','$\sigma = 1.5MPa$','$\sigma = 2MPa$','Interpreter','latex')
title('Square One Element (Force)','Interpreter','latex')

figure()
hold on
plot(u10,d10,'Color',cmp(1,:))
plot(u15,d15,'Color',cmp(2,:))
plot(u20,d20,'Color',cmp(3,:));
legend('$\sigma = 1MPa$','$\sigma = 1.5MPa$','$\sigma = 2MPa$','Interpreter','latex')
title('Square One Element (Damage)','Interpreter','latex')