clc
clear
close all
%% Test to check the performance of the NN
nElem = 10;
th0 = 0.5;
th1 = 2;
xVal = (0:1/nElem:1)';
yS = @(x) th0 + th1*x;

noise = -0.05 + 0.1*rand(length(xVal),1);
y = yS(xVal)+noise; % Fallo aqui

hold on
plot(xVal,y,'o');
plot(xVal,yS(xVal));
hold off

dataset = [xVal, y];

writematrix(dataset, "regression_dataset.csv", 'Delimiter', ",");