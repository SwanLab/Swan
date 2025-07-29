function fittingPhasefield
matType = load('SquareArea.mat');

x = matType.phi;
C11v = squeeze(matType.mat(1,1,1,1,:));
C12v = squeeze(matType.mat(1,1,2,2,:));
C33v = squeeze(matType.mat(1,2,1,2,:));

num = @(p) (p(19).*x.^9 + p(17).*x.^8 + p(15).*x.^7 + p(13).*x.^6 ...
    + p(11).*x.^5 + p(9).*x.^4 + p(7).*x.^3 + p(5).*x.^2 ...
    + p(3).*x.^1  + p(1));
den = @(p) (p(20).*x.^9 + p(18).*x.^8 + p(16).*x.^7 + p(14).*x.^6 ...
    + p(12).*x.^5 + p(10).*x.^4 + p(8).*x.^3 + p(6).*x.^2 ...
    + p(4).*x.^1  + p(2));
C11fit = @(p) num(p)./den(p);

num = @(p) (p(39).*x.^9 + p(37).*x.^8 + p(35).*x.^7 + p(33).*x.^6 ...
    + p(31).*x.^5 + p(29).*x.^4 + p(27).*x.^3 + p(25).*x.^2 ...
    + p(23).*x.^1  + p(21));
den = @(p) (p(40).*x.^9 + p(38).*x.^8 + p(36).*x.^7 + p(34).*x.^6 ...
    + p(32).*x.^5 + p(30).*x.^4 + p(28).*x.^3 + p(26).*x.^2 ...
    + p(24).*x.^1  + p(22));
C12fit = @(p) num(p)./den(p);

num = @(q) (q(19).*x.^9 + q(17).*x.^8 + q(15).*x.^7 + q(13).*x.^6 ...
    + q(11).*x.^5 + q(9).*x.^4 + q(7).*x.^3 + q(5).*x.^2 ...
    + q(3).*x.^1  + q(1));
den = @(q) (q(20).*x.^9 + q(18).*x.^8 + q(16).*x.^7 + q(14).*x.^6 ...
    + q(12).*x.^5 + q(10).*x.^4 + q(8).*x.^3 + q(6).*x.^2 ...
    + q(4).*x.^1  + q(2));
C33fit = @(q) num(q)./den(q);






%% C11 - C12
A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];
nonlcon = @(p) nonLinearConC11(p,C11v,C12v);

objective = @(p) sum(sqrt(((C11fit(p)'-C11v)./C11v).^2) + sqrt(((C12fit(p)'-C12v)./C12v).^2));
options = optimoptions(@fmincon,'StepTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'MaxFunctionEvaluations',10000);

objResCon=100;
for i=1:50
    p0 = rand(1,40);
    [popt,fval,exitflag,output,lambda,grad,hessian] = fmincon(objective,p0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    if fval<objResCon
        objResCon = fval;
        pResCon = popt;
    end
end
disp("Initial objective: " + num2str(objective(p0)));
disp("Final objective: "   + num2str(objective(pResCon)));
plot(x,C11v,'ro')
hold on
plot(x,C11fit(pResCon),'-')
legend('measured','optimal')

figure()
plot(x,C12v,'ro')
hold on
plot(x,C12fit(pResCon),'-')

%% C33 
A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];
nonlcon = @(q) nonLinearConC33(q,C33v);

objective = @(q) sum(sqrt(((C33fit(q)'-C33v)./C33v).^2));
options = optimoptions(@fmincon,'StepTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'MaxFunctionEvaluations',10000);

objResCon=100;
for i=1:1000
    q0 = rand(1,20);
    [qopt,fval,exitflag,output,lambda,grad,hessian] = fmincon(objective,q0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    if fval<objResCon
        objResCon = fval;
        qResCon = qopt;
    end
end
disp("Initial objective: " + num2str(objective(q0)));
disp("Final objective: "   + num2str(objective(qResCon)));
plot(x,C33v,'ro')
hold on
plot(x,C33fit(qResCon),'-')
legend('measured','optimal')

%% Save data
syms x
p = pResCon;
num = p(19).*x.^9 + p(17).*x.^8 + p(15).*x.^7 + p(13).*x.^6 ...
    + p(11).*x.^5 + p(9).*x.^4 + p(7).*x.^3 + p(5).*x.^2 ...
    + p(3).*x.^1  + p(1);
den = p(20).*x.^9 + p(18).*x.^8 + p(16).*x.^7 + p(14).*x.^6 ...
    + p(12).*x.^5 + p(10).*x.^4 + p(8).*x.^3 + p(6).*x.^2 ...
    + p(4).*x.^1  + p(2);
C11fun = num./den;

num = p(39).*x.^9 + p(37).*x.^8 + p(35).*x.^7 + p(33).*x.^6 ...
    + p(31).*x.^5 + p(29).*x.^4 + p(27).*x.^3 + p(25).*x.^2 ...
    + p(23).*x.^1  + p(21);
den = p(40).*x.^9 + p(38).*x.^8 + p(36).*x.^7 + p(34).*x.^6 ...
    + p(32).*x.^5 + p(30).*x.^4 + p(28).*x.^3 + p(26).*x.^2 ...
    + p(24).*x.^1  + p(22);
C12fun = num./den;

q = qResCon;
num = q(19).*x.^9 + q(17).*x.^8 + q(15).*x.^7 + q(13).*x.^6 ...
    + q(11).*x.^5 + q(9).*x.^4 + q(7).*x.^3 + q(5).*x.^2 ...
    + q(3).*x.^1  + q(1);
den = q(20).*x.^9 + q(18).*x.^8 + q(16).*x.^7 + q(14).*x.^6 ...
    + q(12).*x.^5 + q(10).*x.^4 + q(8).*x.^3 + q(6).*x.^2 ...
    + q(4).*x.^1  + q(2);
C33fun = num./den;

mat = matType.mat;
phi = matType.phi;
degradation = matType.degradation;
degradation.fun{1,1,1,1} = matlabFunction(C11fun);
degradation.fun{1,1,2,2} = matlabFunction(C12fun);
degradation.fun{2,2,1,1} = matlabFunction(C12fun);
degradation.fun{2,2,2,2} = matlabFunction(C11fun);
degradation.fun{1,2,1,2} = matlabFunction(C33fun);
degradation.fun{2,1,2,1} = matlabFunction(C33fun);

degradation.dfun{1,1,1,1} = matlabFunction(diff(C11fun));
degradation.dfun{1,1,2,2} = matlabFunction(diff(C12fun));
degradation.dfun{2,2,1,1} = matlabFunction(diff(C12fun));
degradation.dfun{2,2,2,2} = matlabFunction(diff(C11fun));
degradation.dfun{1,2,1,2} = matlabFunction(diff(C33fun));
degradation.dfun{2,1,2,1} = matlabFunction(diff(C33fun));

degradation.ddfun{1,1,1,1} = matlabFunction(diff(diff(C11fun)));
degradation.ddfun{1,1,2,2} = matlabFunction(diff(diff(C12fun)));
degradation.ddfun{2,2,1,1} = matlabFunction(diff(diff(C12fun)));
degradation.ddfun{2,2,2,2} = matlabFunction(diff(diff(C11fun)));
degradation.ddfun{1,2,1,2} = matlabFunction(diff(diff(C33fun)));
degradation.ddfun{2,1,2,1} = matlabFunction(diff(diff(C33fun)));
save('NewSet1s','mat','phi','degradation');
end



function C = createC(x,fun,phi)
p11 = x(1:2:20);
q11 = x(2:2:20);
p12 = x(21:2:40);
q12 = x(22:2:40);

C11  = fun(p11,q11,phi);
C12  = fun(p12,q12,phi);
C    = [C11 C12 0; C12 C11 0; 0 0 1];
end

%% Function

function [J,dJ] = costFunction()


end

function [c,ceq] = nonLinearConC11(x,C11v,C12v)

%E=210; nu=0.3; C=E/((1+nu)*(1-nu));
E = 210;
Gc=5e-3; l0=0.1; sigCrit=1.5; cw=8/3;
K = 2*(Gc/(cw*l0))*E;

p11 = x(1:2:20);
q11 = x(2:2:20);
p12 = x(21:2:40);
q12 = x(22:2:40);


C  = @(phi) createC(x,rational)
dC = @(phi) createC(x,rationalDeriv)

sig = [sigCrit 0 0];

dE = @(phi) sig*inv(C(phi))*dC(phi)*inv(C(phi))*sig;

C11 = C(1,1);
C12 = C(1,2);

ceq(1) = C11(0) - C11v(1);
ceq(2) = C12(0) - C12v(1);
ceq(3) = C11(1) - 0;
ceq(4) = C12(1) - 0;
ceq(5) = dE(0) - K;



c = [];
% ceq(1) = x(1)/x(2) - C11v(1) ;
% ceq(2) = x(21)/x(22) - C12v(1);
% ceq(3) = (x(1)+x(3)+x(5)+x(7)+x(9)+x(11)+x(13)+x(15)+x(17)+x(19))/...
%          (x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)+x(20));
% ceq(4) = (x(21)+x(23)+x(25)+x(27)+x(29)+x(31)+x(33)+x(35)+x(37)+x(39))/...
%          (x(22)+x(24)+x(26)+x(28)+x(30)+x(32)+x(34)+x(36)+x(38)+x(40));



%ceq(5) = 2*(nu^2)*((p(23)*p(22)-p(21)*p(24))/(p(22)^2)) ...
%         -(1+nu^2)*((p(3)*p(2)-p(1)*p(4))/(p(2)^2)) ...
%          -2*(Gc/(cw*l0))*((E^2)/C)*(1/sigCrit)^2;
% 
% ceq(5) = 2*(nu)*((x(23)*x(22)-x(21)*x(24))/(x(22)^2)) ...
%          -(1+nu^2)*((x(3)*x(2)-x(1)*x(4))/(x(2)^2)) ...
%          -2*(Gc/(cw*l0))*((E^2))*(1/sigCrit)^2;

end



function [c,ceq] = nonLinearConC33(q,C33)
c = [];
ceq(1) = q(1)/q(2) - C33(1) ;
ceq(2) = (q(1)+q(3)+q(5)+q(7)+q(9)+q(11)+q(13)+q(15)+q(17)+q(19))/...
         (q(2)+q(4)+q(6)+q(8)+q(10)+q(12)+q(14)+q(16)+q(18)+q(20));

% E=210; nu=0.3; C=E/((1+nu)*(1-nu));
% Gc=5e-3; l0=0.1; sigCrit=1.5; cw=8/3;
end

% 2*(nu^2)*degradation.dfun{1,1,2,2}(0) - (1+nu^2)*degradation.dfun{1,1,1,1}(0) ...
%          -2*(Gc/(cw*l0))*((E^2)/C)*(1/sigCrit)^2;




function y = polyEval(p,phi)
    y = 0;
    n = length(p);
    for i = 1:n
        y = y + p(i).* phi.^i;
    end
end

function y = rational(p,q,phi)
    y = polyEval(p,phi)./polyEval(q,phi);
end


function y = polyDeriv(p,phi)
    y = 0;
    n = length(p);
    for i = 1:n
        if i > 1
            y = y + (i) * p(i) .* phi.^(i-1);
        end
    end
end

function y = rationalDeriv(p,q,phi)
    N  = polyEval(p,phi);
    D  = polyEval(q,phi);
    Np = polyDeriv(p,phi);
    Dp = polyDeriv(q,phi);

    y = (Np .* D - N .* Dp) ./ (D.^2);
end
