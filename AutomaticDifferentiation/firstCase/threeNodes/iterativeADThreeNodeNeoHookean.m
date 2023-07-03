function [val, grad] = iterativeADThreeNodeNeoHookean(u0)

C1 = 1; %material constant
D1 = 1; %material constant
D = [1,-1,0; -1,2,-1; 0,-1,1];

[~,numElem] = size(u0);

G = zeros(1,numElem);
G(1,numElem) = 1 * 10^(0);

u1 = ValGradForward(u0(1),[1 0 0]);
u2 = ValGradForward(u0(2),[0 1 0]);
u3 = ValGradForward(u0(3),[0 0 1]);

F(1,1) = D(1,1) * u1; F(1,2) = D(1,2) * u2; F(1,3) = D(1,3) * u3;
F(2,1) = D(1,2) * u1; F(2,2) = D(2,2) * u2; F(2,3) = D(2,3) * u3;
F(3,1) = D(1,3) * u1; F(3,2) = D(2,3) * u2; F(3,3) = D(3,3) * u3;

I1 = F(1,1)^2 + F(2,2)^2+ F(3,3)^2;

% J = F(1,1) * ( F(2,2) * F(3,3) - F(2,3) * F(3,2) ) - F(1,2) * ( F(2,1) * F(3,3) - F(2,3) * F(3,1) ) + F(1,3) * ( F(2,1) * F(3,2) - F(2,2) * F(3,1) );

J = F(1,1) * F(2,2) * F(3,3);

fun1 = C1 * ( I1 - 3 );

fun2 = 2 * C1 * log(J);

fun3 = D1 * ( J - 1 )^2;

fun4 = G(1) * u1 + G(2) * u2 + G(3) * u3;

fun = fun1 - fun2 + fun3 - fun4;

AD = fun.double;
val = AD(:,1);
grad = AD(:,2:end);