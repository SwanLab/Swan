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

%% AUXETIC
DOFs = [200016 399376 747636 1723200 2488840];

iter_ILU = [3629 5437 6793 8728 13197];
iter_AMG = [177 209 214 208 216];
iter_EIFEM_C = [1571 2374 2954 3219 4018];
iter_EIFEM_D = [75 75 76 76 77];

figure;
plot(DOFs,iter_ILU,'-s','LineWidth',2,'MarkerSize',8); hold on;
plot(DOFs,iter_AMG,'-d','LineWidth',2,'MarkerSize',8);
plot(DOFs,iter_EIFEM_C,'-^','LineWidth',2,'MarkerSize',8);
plot(DOFs,iter_EIFEM_D,'-o','LineWidth',2,'MarkerSize',8);

legend('ILU','AMG','EIFEM Continuous','EIFEM Discontinuous P2');
title('Iter. vs DOFs');
xlabel('DOFs');
ylabel('Iterations');
% set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
