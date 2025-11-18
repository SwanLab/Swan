function fittingLessDenominator(degradationData,matInfo,coeffNumber,outputName)
    matType = load(degradationData);
    phiData = matType.phi;
    g11data = squeeze(matType.mat(1,1,1,1,:));
    g12data = squeeze(matType.mat(1,1,2,2,:));
    g33data = squeeze(matType.mat(1,2,1,2,:));
    gData = [g11data';g12data';g33data'];
    
    A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];
    nonlcon = @(coeff) nonLinearCon(coeff,gData,coeffNumber,matInfo);
    objective = @(p) objectiveFun(p,coeffNumber,phiData,gData);
    options = optimoptions(@fmincon, ...
                       'StepTolerance',1e-10, ...
                       'OptimalityTolerance',1e-10, ...
                       'MaxFunctionEvaluations',10000, ...
                       'Display','none');

    % load('Degradation/L0_variation/DegSqr15lHS.mat')
    % [~, ceq] = nonlcon(coeffBestResult);
    % disp(['MinObjective: ', num2str(objBestResult, '%.2e'), ' | MinCeq: [', num2str(abs(ceq), '%.2e '), ']']);
    
    nTries = 250;
    coeff = zeros(nTries,60);
    cost  = zeros(nTries,1);

    parfor_progress(nTries);
    parfor i=1:nTries
        coeff0 = rand(1,60);
        [coeff(i,:),cost(i),~,~,~,~,~] = fmincon(objective,coeff0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        parfor_progress;
    end
    parfor_progress(0);

    [costOpt,tryOpt] = min(cost);
    coeffOpt = coeff(tryOpt,:);

    [~, ceq] = nonlcon(coeffOpt);
    disp(['MinObjective: ', num2str(costOpt, '%.2e'), ' | MinCeq: [', num2str(abs(ceq), '%.2e '), ']']);

    degFuns = generateConstitutiveTensor(coeffOpt,coeffNumber,costOpt,matType,outputName);
    plotResults(objective,phiData,gData,coeffOpt,degFuns)
end


%% Main functions

function [c,ceq] = nonLinearCon(coeff,gData,nCoeff,mat)
    gFun  = @(phi) degradationTensor(coeff,nCoeff,phi);
    dgFun = @(phi) degradationTensorDerivative(coeff,nCoeff,phi);
   
    c = [];
    g0 = cell2mat(gFun(0));
    ceq(1) = g0(1,1) - gData(1,1);
    ceq(2) = g0(1,2) - gData(2,1);
    ceq(3) = g0(3,3) - gData(3,1);
    g1 = cell2mat(gFun(1));
    ceq(4) = g1(1,1);
    ceq(5) = g1(1,2);
    ceq(6) = g1(3,3);

    cw=8/3;
    sigCrit=[mat.sigma 0 0];
    
    dg0 = cell2mat(dgFun(0));
    ceq(7) = (1/2)*sigCrit*((g0\dg0)/g0)*sigCrit' + mat.E*mat.Gc/(cw*mat.l0);
    %dC1 = cell2mat(dCfun(1));
    % ceq(8)  = dC1(1,1);
    % ceq(9)  = dC1(1,2);
    % ceq(10) = dC1(3,3);
end

function objFun = objectiveFun(coeff,nCoeffs,phiData,Cdata)
    C  = degradationTensor(coeff,nCoeffs,phiData);
    idx = [1 1; 1 2; 3 3];
    objFun = 0;
    for i=1:3
        objFun = objFun + sum(((C{idx(i,1),idx(i,2)}- Cdata(i,:))./Cdata(i,:)).^2);
    end
end

%% Definition of functions

function g = degradationTensor(coeff,nCoeffs,phi)
    x = cell(1,length(nCoeffs));
    for i=1:length(nCoeffs)
        x{i} = coeff(sum(nCoeffs(1:i-1))+1:sum(nCoeffs(1:i)));
    end
    g11fun = rationalFun(x{1},x{2},phi);
    g12fun = rationalFun(x{3},x{4},phi);
    g33fun = rationalFun(x{5},x{6},phi); 
    
    g = {g11fun, g12fun,   0;
         g12fun, g11fun,   0;
           0,      0,   g33fun};
end

function fun = rationalFun(p,q,phi)
    num = polyFun(p,phi);
    den = polyFun(q,phi);
    fun = num./den;
end

function fun = polyFun(a,phi)
    n = length(a);
    fun = a(1);
    for i=1:n-1
        fun = fun + a(i+1)*phi.^(i);
    end
end

function dC = degradationTensorDerivative(coeff,nCoeffs,phi)
    x = cell(1,length(nCoeffs));
    for i=1:length(nCoeffs)
        x{i} = coeff(sum(nCoeffs(1:i-1))+1:sum(nCoeffs(1:i)));
    end
    dg11fun = rationalFunDeriv(x{1},x{2},phi);
    dg12fun = rationalFunDeriv(x{3},x{4},phi);
    dg33fun = rationalFunDeriv(x{5},x{6},phi); 
    
    dC = {dg11fun, dg12fun,   0;
          dg12fun, dg11fun,   0;
             0,       0,   dg33fun};
end

function fun = rationalFunDeriv(p,q,phi)
    N  = polyFun(p,phi);
    dN = polyFunDeriv(p,phi);
    D  = polyFun(q,phi);
    dD = polyFunDeriv(q,phi);
    fun = (dN.*D - N.*dD)./(D.^2);
end

function fun = polyFunDeriv(a,phi)
    n = length(a);
    fun = a(2);
    for i=2:n-1
        fun = fun + (i)*a(i+1)*phi.^(i-1);
    end
end

%% Others

function plotResults(objective,phiData,Cdata,coeffRes,degFuns)
    close all
    disp("Final objective: "   + num2str(objective(coeffRes)));

    figure()
    plot(phiData,Cdata(1,:),'ro')
    hold on
    fplot(degFuns.fun{1,1,1,1},[0 1],'-')
    title('C11')

    figure()
    plot(phiData,Cdata(2,:),'ro')
    hold on
    fplot(degFuns.fun{1,1,2,2},[0 1],'-')
    title('C12')

    figure()
    plot(phiData,Cdata(3,:),'ro')
    hold on
    fplot(degFuns.fun{1,2,1,2},[0 1],'-')
    legend('measured','optimal')
    title('C33')
end

function degradation = generateConstitutiveTensor(coeffOpt,nCoeffs,costOpt, matType, name)
    [C11, C12, C33] = recoverTensorComponents(coeffOpt,nCoeffs);
    mat = matType.mat;
    phi = matType.phi;
    degradation = matType.degradation;

    comps = {C11, C12, C33};
    idx = {{1,1,1,1}, {1,1,2,2}, {2,2,1,1}, {2,2,2,2}, {1,2,1,2}, {2,1,2,1}};
    map = [1 2 2 1 3 3];

    for i = 1:length(idx)
        degradation.fun{idx{i}{:}} = matlabFunction(comps{map(i)});
        degradation.dfun{idx{i}{:}} = matlabFunction(diff(comps{map(i)}));
        degradation.ddfun{idx{i}{:}} = matlabFunction(diff(diff(comps{map(i)})));
    end

    save(name, 'mat', 'phi', 'degradation','coeffOpt','costOpt');
end


function [C11,C12,C33] = recoverTensorComponents(coeff,nCoeffs)
    x = cell(1,length(nCoeffs));
    for i=1:length(nCoeffs)
        x{i} = coeff(sum(nCoeffs(1:i-1))+1:sum(nCoeffs(1:i)));
    end
    syms phi
    C11 = rationalFun(x{1},x{2},phi);
    C12 = rationalFun(x{3},x{4},phi);
    C33 = rationalFun(x{5},x{6},phi); 
end