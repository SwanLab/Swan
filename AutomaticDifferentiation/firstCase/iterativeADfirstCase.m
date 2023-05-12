function [val, grad] = iterativeADfirstCase(u)

D = [1,-1,1;-1,1,1;1,1,1];

u1 = ValGradForward(u(1),[1 0 0]);
u2 = ValGradForward(u(2),[0 1 0]);
u3 = ValGradForward(u(3),[0 0 1]);

W = 0.5 * ( ( D(1,1)^2 + D(2,1)^2 + D(3,1)^2 ) * u1^2 + ( D(1,2)^2 + D(2,2)^2 + D(3,2)^2 ) * u2^2 + ( D(1,3)^2 + D(2,3)^2 + D(3,3)^2 ) * u3^2 + ( D(1,1) * D(1,2) + D(2,1) * D(2,2) + D(3,1) * D(3,2) ) * u1 * u2 + ( D(1,1) * D(1,3) + D(2,1) * D(2,3) + D(3,1) * D(3,3) ) * u1 * u3 + ( D(1,2) * D(1,3) + D(2,2) * D(2,3) + D(3,2) * D(3,3) ) * u2 * u3);

AD = W.double;
val = AD(:,1);
grad = AD(:,2:end);