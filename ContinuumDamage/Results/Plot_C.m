close all;

load('Csec.mat');
load('Ctan.mat');
load('DispVal.mat');

figure();
hold on
plot(DispVal,Ctan,'DisplayName','Ctan');
plot(DispVal,Csec,'DisplayName','Csec');

hold off
legend; 
title('Material - Displacement');