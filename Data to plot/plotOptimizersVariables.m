%% OPTIMIZERS VARIABLES COMPARISON 40x40 MESH %%
close all
load("AugmentedLagrCant04.mat");
cost.Lag  = c;
it.Lag    = 1:length(c);
% merit.Lag = m;
const.Lag = h;
% tau.Lag   = t;
grad.Lag  = g;
dual.Lag  = d;

load("NullSpaceCant04.mat");
cost.Null  = c(1:200);
it.Null    = 1:200;
% merit.Null = m;
const.Null = h(1:200);
% tau.Null   = t;
grad.Null  = g(1:200);
dual.Null  = d(1:200);
var.Null   = v(1:200);

% load("NullSpaceNoConvergence.mat");
% cost.Null2  = c;
% it.Null2    = 1:length(cost.Null2);
% merit.Null2 = m;
% const.Null2 = h;
% tau.Null2   = t;
% grad.Null2  = g;
% dual.Null2   = d;

load("BisectionVariablesCant04.mat");
cost.Bis  = c;
it.Bis    = 1:length(c);
% merit.Bis = m;
const.Bis = h;
% tau.Bis   = t;
grad.Bis  = g;
dual.Bis  = d;

figure(1)
plot(it.Lag,cost.Lag,"red",it.Null,cost.Null,"blue",it.Bis,cost.Bis,"green")
xlabel('Iterations')
ylabel('Cost')
legend('Augmented Lagrangian', 'Null Space','Null Space', 'Bisection')
figure(2)
plot(it.Lag,const.Lag,"red",it.Null,const.Null,"blue",it.Bis,const.Bis,"green")
xlabel('Iterations')
ylabel('Constraint')
legend('Augmented Lagrangian', 'Null Space','Null Space', 'Bisection')
% legend('Augmented Lagrangian', 'Null Space','Null Space w/o \tau restriction', 'Bisection')
figure(3)
% plot(it.Lag,merit.Lag,"red",it.Null,merit.Null,"blue",it.Bis,merit.Bis,"green")
% title('Merit function value')
% % legend('Augmented Lagrangian', 'Null Space','Null Space w/o \tau restriction', 'Bisection')
% subplot(3,3,4)
% plot(it.Lag,tau.Lag,"red",it.Null,tau.Null,"blue",it.Bis,tau.Bis,"green")
% title('Step length')
% legend('Augmented Lagrangian', 'Null Space','Null Space w/o \tau restriction', 'Bisection')
% subplot(3,3,5)
plot(it.Lag,grad.Lag,"red",it.Null,grad.Null,"blue",it.Bis,grad.Bis,"green")
xlabel('Iterations')
ylabel('Cost gradient norm')
legend('Augmented Lagrangian', 'Null Space','Null Space', 'Bisection')
% legend('Augmented Lagrangian', 'Null Space','Null Space w/o \tau restriction', 'Bisection')
% subplot(3,3,6)
figure(4)
plot(it.Lag,dual.Lag,"red",it.Null,dual.Null,"blue",it.Bis,dual.Bis,"green")
xlabel('Iterations')
ylabel('Dual variable \lambda')
legend('Augmented Lagrangian', 'Null Space','Null Space', 'Bisection')
% legend('Augmented Lagrangian', 'Null Space','Null Space w/o \tau restriction', 'Bisection')


