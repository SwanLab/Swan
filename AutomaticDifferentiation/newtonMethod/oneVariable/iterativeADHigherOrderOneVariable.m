function [val, grad, grad2] = iterativeADHigherOrderOneVariable(u)
x = ValGradForward(ValGradForward(u,1),1);
f = exp(-sqrt(x))*sin(x*log(1+x^2));

firstGrad = f.val.double;
secondGrad = f.grad.double;

val = firstGrad(:,1);
grad = firstGrad(:,2);
grad2 = secondGrad(:,2);