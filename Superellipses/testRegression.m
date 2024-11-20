clc
clear
close all
%% Test to check the performance of the NN
nElem = 10;
th0 = 0.5;
th1 = 2;
xVal = (0:1/nElem:1)';
yS = @(x) th0 + th1*x;

noise = -0.1 + 0.2*rand(length(xVal),1);
y = yS(xVal)+noise; 

hold on
plot(xVal,y,'o');
plot(xVal,yS(xVal));
hold off
grid minor
box on
ylim([0 2.5])
xlabel('x')
ylabel('y')

dataset = [xVal, y];

writematrix(dataset, "noisyRegression.csv", 'Delimiter', ",");