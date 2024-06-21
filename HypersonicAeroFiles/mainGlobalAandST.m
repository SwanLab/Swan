clc
clear
close all

%% Generalization of the code to allow any a-values and any source term
% Definition of advection velocity function handle
s.aType = 'basic'; % 'vortex' 'stagnation'
s.stType = 'step';

s.method = 3; % Only TG2 method used
s.revolutions = 1;
s.timeSteps = 120;
probM1 = UnsteadyConvectionProblem(s);
u  = probM1.compute();
nLast = size(u,2);
probM1.plot(u,1);
probM1.plot(u,25);
probM1.plot(u,50);
probM1.plot(u,100);
probM1.plot(u,120);
