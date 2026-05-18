%% CIRCLE

% scalability
x=[43968, 116864, 262296, 654480];

dat03     = [0.25, 0.80, 2.37, 7.51];
NN03      = [0.30, 0.91, 2.76, 8.99];
Direct03  = [0.10, 0.33, 0.88, 2.51];

dat08     = [0.60, 2.23, 6.40, 20.74];
NN08      = [0.65, 2.77, 8.63, 31.68];
Direct08  = [0.11, 0.34, 0.87, 2.54];

figure
plot(x,dat03,'-+','LineWidth',1.5);
hold on
plot(x,NN03,'-+','LineWidth',1.5);
hold on
plot(x,Direct03,'-+','LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
% set(gca, 'YScale', 'log')
title('Time vs DOFs for r=0.3');
xlabel('DOFs');
ylabel('time [s]');
ylim([0,15]);

figure
plot(x,dat08,'-+','LineWidth',1.5);
hold on
plot(x,NN08,'-+','LineWidth',1.5);
hold on
plot(x,Direct08,'-+','LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
% set(gca, 'YScale', 'log')
title('Time vs DOFs for r=0.8');
xlabel('DOFs');
ylabel('time [s]');
ylim([0,45]);


% Mesh generalization
x= [10108, 24034, 43968, 101828, 183688, 289548];

dat=[0.05, 0.16, 0.33, 1.07, 2.81, 5.11];
NN=[0.05, 0.14, 0.32, 1.10, 2.81, 5.28];
Direct=[0.02, 0.06, 0.17, 0.26, 0.53, 0.91];

figure
plot(x,dat,'-+','LineWidth',1.5);
hold on
plot(x,NN,'-+','LineWidth',1.5);
hold on
plot(x,Direct,'-+','LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
% set(gca, 'YScale', 'log')
title('Time vs DOFs');
xlabel('DOFs');
ylabel('time [s]');
ylim([0,8]);


%% LATTICE
x=[101828,270944,608616,811376];

dat=[1.17,3.94,9.82,16.70];
NN=[1.29,4.19,10.98,17.15];
Direct=[0.28,0.78,2.16,3.44];

figure
plot(x,dat,'-+','LineWidth',1.5);
hold on
plot(x,NN,'-+','LineWidth',1.5);
hold on
plot(x,Direct,'-+','LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
set(gca, 'YScale', 'log')
title('Time vs DOFs');
xlabel('DOFs');
ylabel('time [s]');
ylim([10^-1 10^2]);

% Mesh generalization
x=[10108,24038,43968,101828,183688,289548];

dat=[0.08,0.25,0.65,2.72,7.91,15.98];
NN=[0.06,0.23,0.67,2.39,6.59,13.08];
Direct=[0.03,0.06,0.10,0.26,0.56,0.90];

figure
plot(x,dat,'-+','LineWidth',1.5);
hold on
plot(x,NN,'-+','LineWidth',1.5);
hold on
plot(x,Direct,'-+','LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
set(gca, 'YScale', 'log')
title('Time vs DOFs');
xlabel('DOFs');
ylabel('time [s]');
ylim([10^-1 10^2]);

%% SPHERE

% Scalability
x1=[264432,457284,775872,1029120];
x2=[264432,457284,775872];

dat=[3.15,7.67,13.19,18.63];
NN=[3.97,8.12,14.5,20.42];
Direct=[54.29,292.63,1048.8];

figure
plot(x1,dat,'LineWidth',1.5);
hold on
plot(x1,NN,'LineWidth',1.5);
hold on
plot(x2,Direct,'LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
set(gca, 'YScale', 'log')
title('Time vs DOFs');
xlabel('DOFs');
ylabel('time [s]');

% Mesh generalization
x1=[11952,82092,264432,612972,924336];
x2=[11952,82092,264432,612972];

dat=[0.09,0.7,3.65,11.96,20.07];
NN=[0.08,0.72,4.04,11.46,19.86];
Direct=[0.1341,4.20,60.54,500.55];

figure
plot(x1,dat,'LineWidth',1.5);
hold on
plot(x1,NN,'LineWidth',1.5);
hold on
plot(x2,Direct,'LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
set(gca, 'YScale', 'log')
title('Time vs DOFs');
xlabel('DOFs');
ylabel('time [s]');


%% CUBE

% Scalability
x1=[264432,457284,775872,1029120];
x2=[264432,457284,775872];

dat=[3.74,7.38,16.83,18.54];
NN=[3.81,8.20,14.21,19.96];
Direct=[48.72,450.72,1361];

figure
plot(x1,dat,'LineWidth',1.5);
hold on
plot(x1,NN,'LineWidth',1.5);
hold on
plot(x2,Direct,'LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
set(gca, 'YScale', 'log')
title('Time vs DOFs');
xlabel('DOFs');
ylabel('time [s]');

% Mesh generalization
x1=[11952,82092,264432,612972,924336];
x2=[11952,82092,264432,612972];

dat=[0.12,0.80,5.10,15.89,30.80];
NN=[0.18,0.96,5.28,16.26,29.09];
Direct=[0.12,8.21,45.84,513.31];

figure
plot(x1,dat,'LineWidth',1.5);
hold on
plot(x1,NN,'LineWidth',1.5);
hold on
plot(x2,Direct,'LineWidth',1.5);
legend('Multiscale', 'Multiscale with NN', 'Direct Solver');
set(gca, 'YScale', 'log')
title('Time vs DOFs');
xlabel('DOFs');
ylabel('time [s]');