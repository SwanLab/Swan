clear;
clc;
close all;

prob = TopOptTestTutorialDensityNullSpace();
rho  = prob.designVariable;
mesh = prob.mesh;
pT   = 0.2;

epsilons = [1 2 3 4 5 6];

Results = cell(2*size(epsilons,2),4);


for ii = 1:size(epsilons,2)
    % p = 1 case
    s.mesh            = mesh;
    s.perimeterTarget = pT;
    s.p               = 1;
    s.eps             = epsilons(ii);
    s.gradientTest    = LagrangianFunction.create(mesh,1,'P1');
    constraintP1_apr1 = PerimeterNormPFunctionalTest(s);
    [J,~]             = constraintP1_apr1.computeFunctionAndGradient(rho);
    perOverVolP1_apr1 = (J+1)*pT;
    constraintP1_apr2 = PerimeterNormPFunctional(s);
    [J,~]             = constraintP1_apr2.computeFunctionAndGradient(rho);
    perOverVolP1_apr2 = (J+1)*pT;
    
    % p = "inf" case
    s.mesh              = mesh;
    s.perimeterTarget   = pT;
    s.p                 = 128;
    s.eps               = epsilons(ii);
    s.gradientTest      = LagrangianFunction.create(mesh,1,'P1');
    constraintPInf_apr1 = PerimeterNormPFunctionalTest(s);
    [J,~]               = constraintPInf_apr1.computeFunctionAndGradient(rho);
    perOverVolPInf_apr1 = (J+1)*pT;
    constraintPInf_apr2 = PerimeterNormPFunctional(s);
    [J,~]               = constraintPInf_apr2.computeFunctionAndGradient(rho);
    perOverVolPInf_apr2 = (J+1)*pT;

    Results{2*ii-1,1}   = '1st Approach';
    Results{2*ii-1,2}   = epsilons(ii);
    Results{2*ii-1,3}   = perOverVolP1_apr1;
    Results{2*ii-1,4}   = perOverVolPInf_apr1;
    Results{2*ii,1}     = '2nd Approach';
    Results{2*ii,2}     = epsilons(ii);
    Results{2*ii,3}     = perOverVolP1_apr2;
    Results{2*ii,4}     = perOverVolPInf_apr2;
end