function [val, grad] = iterativeADNewton(u)
x = ValGradForward(u,1);
f = exp(-sqrt(x))*sin(x*log(1+x^2));

AD = f.double;
val = AD(:,1);
grad = AD(:,2:end);