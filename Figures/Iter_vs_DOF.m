%% CIRCLE
pos= [543 498 433 315]

% SCALABILITY
dof = [43968, 116864, 262296, 654480];

DS_03 = [18, 19, 20, 20];
NN_03 = [19, 22, 24, 24];

DS_08 = [50, 57, 57, 57];
NN_08 = [52, 71, 76, 86];

figure
plot(dof,DS_03,'LineWidth',1.5);
hold on
plot(dof,NN_03,'LineWidth',1.5);
legend('Dataset', 'NN');
% set(gca, 'YScale', 'log')
title('Scalability test: iter. vs DOFs for r = 0.3');
xlabel('DOFs');
ylabel('Iterations');
ylim([10,30]);


figure
plot(dof,DS_08,'LineWidth',1.5);
hold on
plot(dof,NN_08,'LineWidth',1.5);
legend('Dataset', 'NN');
% set(gca, 'YScale', 'log')
title('Scalability test: iter. vs DOFs for r = 0.8');
xlabel('DOFs');
ylabel('Iterations');
ylim([45,95]);


% GENERALIZATION
dof = [10108, 24034, 43968, 101828, 183688, 576240];
DS= [13, 18, 24, 34, 44, 55];
NN= [13, 19, 24, 34, 45, 56];

plot(dof,DS,'LineWidth',1.5);
hold on
plot(dof,NN,'LineWidth',1.5);
legend('Dataset', 'NN');
% set(gca, 'YScale', 'log')
title('Iter. vs DOFs');
xlabel('DOFs');
ylabel('Iterations');


