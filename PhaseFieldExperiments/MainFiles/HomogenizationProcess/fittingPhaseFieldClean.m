function mat = fittingPhaseFieldClean()
    matType = load('CircleArea.mat');
    
    phiData = matType.phi;
    C11data = squeeze(matType.mat(1,1,1,1,:));
    C12data = squeeze(matType.mat(1,1,2,2,:));
    C33data = squeeze(matType.mat(1,2,1,2,:));
    Cdata = [C11data';C12data';C33data'];
    
    A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];
    nonlcon = @(coeff) nonLinearCon(coeff,Cdata);
    objective = @(p) objectiveFun(p,phiData,Cdata);
    options = optimoptions(@fmincon,'StepTolerance',1e-10,'OptimalityTolerance',1e-10,...
                           'MaxFunctionEvaluations',10000);
    objBestResult=100;
    for i=1:25
        coeff0 = rand(1,60);
        [coeffOpt,fval,exitflag,output,lambda,grad,hessian] = fmincon(objective,coeff0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        if i==1 || fval<objBestResult
            objBestResult    = fval;
            coeffBestResult  = coeffOpt;
            coeff0BestResult = coeff0;
        end
    end
    plotResults(objective,phiData,Cdata,coeff0BestResult,coeffBestResult)
    generateConstitutiveTensor(coeffBestResult,matType,'CircleAreaDerivative');
end

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

function [c,ceq] = nonLinearCon(coeff,Cdata)
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

    E=210;
    Gc=5e-3; l0=0.1; cw=8/3;
    sigCrit=[1.5 0 0];
    
    dC0 = cell2mat(dCfun(0));
    ceq(7) = sigCrit*inv(C0)*dC0*inv(C0)*sigCrit' + 2*(Gc/(cw*l0))*(E);
end

function objFun = objectiveFun(coeff,phiData,Cdata)
    C  = constitutiveTensor(coeff,phiData);
    idx = [1 1; 1 2; 3 3];
    objFun = 0;
    for i=1:3
        objFun = objFun + sum(((C{idx(i,1),idx(i,2)}- Cdata(i,:))./Cdata(i,:)).^2);
    end
end

function plotResults(objective,phiData,Cdata,coeff0,coeffRes)
    Cfit = constitutiveTensor(coeffRes,phiData);
    close all
    disp("Initial objective: " + num2str(objective(coeff0)));
    disp("Final objective: "   + num2str(objective(coeffRes)));

    figure()
    plot(phiData,Cdata(1,:),'ro')
    hold on
    plot(phiData,Cfit{1,1},'-')
    title('C11')

    figure()
    plot(phiData,Cdata(2,:),'ro')
    hold on
    plot(phiData,Cfit{1,2},'-')
    title('C12')

    figure()
    plot(phiData,Cdata(3,:),'ro')
    hold on
    plot(phiData,Cfit{3,3},'-')
    legend('measured','optimal')
    title('C33')
end

function generateConstitutiveTensor(coeff,matType,name)
    [C11,C12,C33] = recoverTensorComponents(coeff);
    mat = matType.mat;
    phi = matType.phi;
    degradation = matType.degradation;
    degradation.fun{1,1,1,1} = matlabFunction(C11);
    degradation.fun{1,1,2,2} = matlabFunction(C12);
    degradation.fun{2,2,1,1} = matlabFunction(C12);
    degradation.fun{2,2,2,2} = matlabFunction(C11);
    degradation.fun{1,2,1,2} = matlabFunction(C33);
    degradation.fun{2,1,2,1} = matlabFunction(C33);

    degradation.dfun{1,1,1,1} = matlabFunction(diff(C11));
    degradation.dfun{1,1,2,2} = matlabFunction(diff(C12));
    degradation.dfun{2,2,1,1} = matlabFunction(diff(C12));
    degradation.dfun{2,2,2,2} = matlabFunction(diff(C11));
    degradation.dfun{1,2,1,2} = matlabFunction(diff(C33));
    degradation.dfun{2,1,2,1} = matlabFunction(diff(C33));

    degradation.ddfun{1,1,1,1} = matlabFunction(diff(diff(C11)));
    degradation.ddfun{1,1,2,2} = matlabFunction(diff(diff(C12)));
    degradation.ddfun{2,2,1,1} = matlabFunction(diff(diff(C12)));
    degradation.ddfun{2,2,2,2} = matlabFunction(diff(diff(C11)));
    degradation.ddfun{1,2,1,2} = matlabFunction(diff(diff(C33)));
    degradation.ddfun{2,1,2,1} = matlabFunction(diff(diff(C33)));
    save(name,'mat','phi','degradation');
end

function [C11,C12,C33] = recoverTensorComponents(coeff)
    a=coeff(1:20); b=coeff(21:40); c=coeff(41:60);
    syms phi
    C11 = rationalFun(a,phi);
    C12 = rationalFun(b,phi);
    C33 = rationalFun(c,phi); 
end