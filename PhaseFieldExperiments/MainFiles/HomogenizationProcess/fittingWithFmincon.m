function fittingWithFmincon()
    matType = load('/home/gerard/Documents/GitHub/Swan/src/Problems/Damage/Models/PhaseField/PFVademecum/Degradation/SquareArea.mat');

    phiData = matType.phi;
    C11data = squeeze(matType.mat(1,1,1,1,:));
    C12data = squeeze(matType.mat(1,1,2,2,:));
    C33data = squeeze(matType.mat(1,2,1,2,:));
    Cdata = [C11data';C12data';C33data'];
    
    A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];

    [sigma,E,Gc,l0] = computeProperties();
    nonlcon = @(coeff) nonLinearCon(coeff,Cdata,sigma,E,Gc,l0);
    objective = @(p) objectiveFun(p,phiData,Cdata);
    options = optimoptions(@fmincon, ...
                       'StepTolerance',1e-10, ...
                       'OptimalityTolerance',1e-10, ...
                       'MaxFunctionEvaluations',10000, ...
                       'Display','none');

    objBestResult=10000;
    coeff0BestResult = rand(1,60);

    load('Degradation/L0_variation/DegSqr15lHS.mat')
    [~, ceq] = nonlcon(coeffBestResult);
    disp(['MinObjective: ', num2str(objBestResult, '%.2e'), ' | MinCeq: [', num2str(abs(ceq), '%.2e '), ']']);

    for i=1:250
        coeff0 = rand(1,60);
        [coeffOpt,fval,exitflag,output,lambda,grad,hessian] = fmincon(objective,coeff0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        [~, ceq] = nonlcon(coeffOpt);
        disp(['Step: ',num2str(i),' | Objective: ', num2str(fval, '%.2e'), ' | ceq: [', num2str(abs(ceq), '%.2e '), ']']);        

        if fval<objBestResult
            objBestResult    = fval;
            coeffBestResult  = coeffOpt;
            coeff0BestResult = coeff0;
            disp(['NEW VALUE OBTAINED']);
        end
    end
    [~, ceq] = nonlcon(coeffOpt);
    disp(['MinObjective: ', num2str(objBestResult, '%.2e'), ' | MinCeq: [', num2str(abs(ceq), '%.2e '), ']']);

    degFuns = generateConstitutiveTensor(coeffBestResult,objBestResult,matType,'DegSqr15lHS');
    plotResults(objective,phiData,Cdata,coeff0BestResult,coeffBestResult,degFuns)
end


%% Input

function [sigma,E,Gc,l0] = computeProperties()
    E = 210;
    nu = 0.3;
    Gc0 = 5e-3;
    sigma = 1.5;
    % l0 constant
        % sigmach = 1.5;
        % Gc = Gc0*(sigma^2/sigmach^2);
        % l0 = 0.1;
    
    % l0 variable
        Gc = Gc0; 
        lch = (2*E*Gc)/(sigma^2);

        C11 = E/((1+nu)*(1-nu));
        k  = E./(2.*(1-nu));
        mu = E./(2.*(1+nu));
        slope = (- k - k^2/mu - mu - (2*mu^2 + k)/k)/C11;
        lhs = -2*(3/8)*(Gc/slope)*(E/sigma^2);

        l0 = lhs;
end

%% Functions

function C = constitutiveTensor(coeff,phi)
    a=coeff(1:20); b=coeff(21:40); c=coeff(41:60);
    C11fun = rationalFun(a,phi);
    C12fun = rationalFun(b,phi);
    C33fun = rationalFun(c,phi); 
    
    C = {C11fun, C12fun,   0;
         C12fun, C11fun,   0;
           0,      0,   C33fun};
end

function fun = rationalFun(a,phi)
    p = a(1:2:end);
    q = a(2:2:end);
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

function dC = constitutiveTensorDerivative(coeff,phi)
    a=coeff(1:20); b=coeff(21:40); c=coeff(41:60);
    dC11fun = rationalFunDeriv(a,phi);
    dC12fun = rationalFunDeriv(b,phi);
    dC33fun = rationalFunDeriv(c,phi); 
    
    dC = {dC11fun, dC12fun,   0;
          dC12fun, dC11fun,   0;
             0,       0,   dC33fun};
end

function fun = rationalFunDeriv(a,phi)
    p = a(1:2:end);
    q = a(2:2:end);
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

function [c,ceq] = nonLinearCon(coeff,Cdata,sigma,E,Gc,l0)
    dCfun = @(phi) constitutiveTensorDerivative(coeff,phi);
    Cfun  = @(phi) constitutiveTensor(coeff,phi);
   
    c = [];
    C0 = cell2mat(Cfun(0));
    ceq(1) = C0(1,1) - Cdata(1,1);
    ceq(2) = C0(1,2) - Cdata(2,1);
    ceq(3) = C0(3,3) - Cdata(3,1);
    C1 = cell2mat(Cfun(1));
    ceq(4) = C1(1,1);
    ceq(5) = C1(1,2);
    ceq(6) = C1(3,3);

    cw=8/3;
    sigCrit=[sigma 0 0];
    
    dC0 = cell2mat(dCfun(0));
    ceq(7) = (1/2)*sigCrit*inv(C0)*dC0*inv(C0)*sigCrit' + E*Gc/(cw*l0);
    dC1 = cell2mat(dCfun(1));
    % ceq(8)  = dC1(1,1);
    % ceq(9)  = dC1(1,2);
    % ceq(10) = dC1(3,3);
end

function objFun = objectiveFun(coeff,phiData,Cdata)
    C  = constitutiveTensor(coeff,phiData);
    idx = [1 1; 1 2; 3 3];
    objFun = 0;
    for i=1:3
        objFun = objFun + sum(((C{idx(i,1),idx(i,2)}- Cdata(i,:))./Cdata(i,:)).^2);
    end
end

function plotResults(objective,phiData,Cdata,coeff0,coeffRes,degFuns)
    close all
    disp("Initial objective: " + num2str(objective(coeff0)));
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

function degradation = generateConstitutiveTensor(coeffBestResult, objBestResult, matType, name)
    [C11, C12, C33] = recoverTensorComponents(coeffBestResult);
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

    save(name, 'mat', 'phi', 'degradation','coeffBestResult','objBestResult');
end


function [C11,C12,C33] = recoverTensorComponents(coeff)
    a=coeff(1:20); b=coeff(21:40); c=coeff(41:60);
    syms phi
    C11 = rationalFun(a,phi);
    C12 = rationalFun(b,phi);
    C33 = rationalFun(c,phi); 
end