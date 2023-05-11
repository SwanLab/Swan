function [val, grad] = iterativeAD(coord)

x = ValGradForward(coord(1),[1 0]);
y = ValGradForward(coord(2),[0 1]);

f = (x - 1)^2 + (y - 1)^2;

AD = f.double;
val = AD(:,1);
grad = AD(:,2:end);