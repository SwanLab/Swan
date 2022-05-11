%% OPTIMIZERS VARIABLES COMPARISON 40x40 MESH %%

load("AugLagrangianVariables.mat");
cost.Lag  = c;
it.Lag    = 1:length(c);
merit.Lag = m;
const.Lag = h;
tau.Lag   = t;
grad.Lag  = g;
dual.Lag  = d;

load("NullSpaceVariables.mat");
cost.Null  = c;
it.Null    = 1:length(cost.Null);
merit.Null = m;
const.Null = h;
tau.Null   = t;
grad.Null  = g;
dual.Null   = d;

load("BisectionVariables.mat");
cost.Bis  = c;
it.Bis    = 1:length(c);
merit.Bis = m;
const.Bis = h;
tau.Bis   = t;
grad.Bis  = g;
dual.Bis  = d;

figure(1)
hold on
subplot(3,3,1)
plot(it.Lag,cost.Lag,"red",it.Null,cost.Null,"blue",it.Bis,cost.Bis,"green")
title('Cost')
legend('Augmented Lagrangian', 'Null Space', 'Bisection')
subplot(3,3,2)
plot(it.Lag,const.Lag,"red",it.Null,const.Null,"blue",it.Bis,const.Bis,"green")
title('Constraint')
legend('Augmented Lagrangian', 'Null Space', 'Bisection')
subplot(3,3,3)
plot(it.Lag,merit.Lag,"red",it.Null,merit.Null,"blue",it.Bis,merit.Bis,"green")
title('Merit function value')
legend('Augmented Lagrangian', 'Null Space', 'Bisection')
subplot(3,3,4)
plot(it.Lag,tau.Lag,"red",it.Null,tau.Null,"blue",it.Bis,tau.Bis,"green")
title('Step length')
legend('Augmented Lagrangian', 'Null Space', 'Bisection')
subplot(3,3,5)
plot(it.Lag,grad.Lag,"red",it.Null,grad.Null,"blue",it.Bis,grad.Bis,"green")
title('Cost gradient norm')
legend('Augmented Lagrangian', 'Null Space', 'Bisection')
subplot(3,3,6)
plot(it.Lag,dual.Lag,"red",it.Null,dual.Null,"blue",it.Bis,dual.Bis,"green")
title('Dual variable \lambda')
legend('Augmented Lagrangian', 'Null Space', 'Bisection')
hold off

