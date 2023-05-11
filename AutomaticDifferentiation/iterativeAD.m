function [val, grad] = iterativeAD(x)

x1 = ValDerForward(x(1),[1 0 0]);
x2 = ValDerForward(x(2),[0 1 0]);
x3 = ValDerForward(x(3),[0 0 1]);

f = ( x1 * x2 * sin(x3) + exp( x1 * x2 )) / x3;
%f =  sin(x1) * sin(x2) * sin(x3);

AD = f.double;
val = AD(:,1);
grad = AD(:,2:end);