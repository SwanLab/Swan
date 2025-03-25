
clear;
clc;
close all;

prob = TopOptTestTutorialDensityNullSpace();
rho  = prob.designVariable;
mesh = prob.mesh;
pT   = 0.2;

% p = 1 case
s.mesh            = mesh;
s.perimeterTarget = pT;
s.p               = 1;
s.eps             = 2;
s.gradientTest    = LagrangianFunction.create(mesh,1,'P1');
constraintP1      = PerimeterNormPFunctional(s);
[J,~]             = constraintP1.computeFunctionAndGradient(rho);
perOverVolP1      = (J+1)*pT;

% p = "inf" case
s.mesh            = mesh;
s.perimeterTarget = pT;
s.p               = 128;
s.eps             = 2;
s.gradientTest    = LagrangianFunction.create(mesh,1,'P1');
constraintPInf    = PerimeterNormPFunctional(s);
[J,~]             = constraintPInf.computeFunctionAndGradient(rho);
perOverVolPInf    = (J+1)*pT;