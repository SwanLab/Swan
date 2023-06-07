function [val, grad, grad2] = iterativeADNewtonNVariablesNeohookean(u0)

[~,numElem] = size(u0);

C1 = 1;
D1 = 1;
D = eye(numElem);
G = zeros(1,numElem);
G(1,numElem) = 1;

firstGrad = zeros(1,numElem);
firstGrad(1,1) = 1;

secondGrad = zeros(numElem);
secondGrad(1,:) = 1;

u(1) = ValGradForward(ValGradForward(0, firstGrad), secondGrad);

for i = 2:numElem

    firstGrad = zeros(1,numElem);
    firstGrad(1,i) = 1;

    secondGrad = zeros(numElem);
    secondGrad(i,:) = 1;

    u(i) = ValGradForward(ValGradForward(u0(i), firstGrad), secondGrad);
end

% F(1,1) = D(1,1) * u1; F(1,2) = D(1,2) * u2; F(1,3) = D(1,3) * u3;
% F(2,1) = D(1,2) * u1; F(2,2) = D(2,2) * u2; F(2,3) = D(2,3) * u3;
% F(3,1) = D(1,3) * u1; F(3,2) = D(2,3) * u2; F(3,3) = D(3,3) * u3;

for i = 1:numElem

    for j = 1:numElem

        F(i,j) = D(i,j) * u(j);

    end

end

% I1 = F(1,1) + F(2,2) + F(3,3);

I1 = F(1,1);

for i = 2:numElem
    I1 = I1 + F(i,i);
end

J = F(1,1) * ( F(2,2) * F(3,3) - F(2,3) * F(3,2) ) - F(1,2) * ( F(2,1) * F(3,3) - F(2,3) * F(3,1) ) + F(1,3) * ( F(2,1) * F(3,2) - F(2,2) * F(3,1) );

%J = det(F);

fun1 = C1 * ( I1 - 2 );

fun2 = 2 * C1 * log(J);

fun2.val.val(isinf(fun2.val.val)) = -10^(-3);
fun2.val.grad(isinf(fun2.val.grad)) = -10^(-3);
fun2.grad.val(isinf(fun2.grad.val)) = -10^(-3);
fun2.grad.grad(isinf(fun2.grad.grad)) = -10^(-3);

fun2.val.val(isnan(fun2.val.val)) = 0;
fun2.val.grad(isnan(fun2.val.grad)) = 0;
fun2.grad.val(isnan(fun2.grad.val)) = 0;
fun2.grad.grad(isnan(fun2.grad.grad)) = 0;

fun3 = D1 * ( J - 1 )^2;
 
% fun4 = G(1) * u1 + G(2) * u2 + G(3) * u3;

force = G(1) * u(1);

for i = 2:numElem
    force = force + G(i) * u(i);
end

fun = fun1 - fun2 + fun3 - force;

firstGrad = fun.val.double;
secondGrad = fun.grad.double;

val = firstGrad(:,1);
grad = firstGrad(:,2:end);
grad2 = secondGrad(1:end,(1+numElem):end);