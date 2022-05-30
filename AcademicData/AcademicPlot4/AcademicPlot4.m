%% ACADEMIC PLOT TEST 4 %%
close all

load("NullSpaceVariablesAcademicProva.mat");
cost.Null  = c;
it.Null    = 1:length(cost.Null);
% merit.Null = m;
const.Null = h;
% tau.Null   = t;
grad.Null  = g;
dual.Null  = d;
var.Null   = v;
k = var.Null;

load("AugmentedLagrAcademic4.mat");
k2 = v;

load("fminconSQPacademic4.mat");
k3 = v;

load("fminconIPOPTacademic4.mat");
k4 = v;

load("AugmentedLagrAcademic4-2.mat");
k5 = v;

% load()

% c1 = @(x) 1/x(1) - x(2);
% c2 = @(x) x(1)+x(2)-3;

v = -3:0.01:7;  % plotting range from -5 to 5
[x,y] = meshgrid(v);  % get 2-D mesh for x and y
cond1 = -x.^2 + y < 0;  % check conditions for these values
cond2 = -x-y-2 < 0;
cond1 = double(cond1);  % convert to double for plotting
cond2 = double(cond2);
cond1(cond1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
cond2(cond2 == 0) = NaN;
cond = cond1.*cond2;  % multiply the two condaces to keep only the common points
z = 10*ones(length(v(:,2)));
figure()
pcolor(x,y,cond);
shading interp
hold on
plot(k(1,:),k(2,:),'r','LineWidth',2)
plot(k2(1,:),k2(2,:),'b','LineWidth',2)
plot(k3(1,:),k3(2,:),'g','LineWidth',2)
plot(k4(1,:),k4(2,:),'m','LineWidth',2)
plot(k4(1,1),k4(2,1),'o','Color','k','LineWidth',2)
plot(k4(1,end),k4(2,end),'p','Color','k','LineWidth',2)
% plot(k5(1,:),k5(2,:),'k','LineWidth',1)
grid on
legend('','Null Space','Augmented Lagrangian','fmincon-SQP','fmincon-IPOPT','Initial point','Optimum')
title('Academic test 4')
xlabel('x_1')
ylabel('x_2')
% view(0,90)    % change to top view
hold off