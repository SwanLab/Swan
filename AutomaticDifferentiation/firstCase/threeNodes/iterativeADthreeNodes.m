function [val, grad] = iterativeADthreeNodes(u)

G = [0; 0; 1];
G = transpose(G);

u1 = ValGradForward(0, [1 0 0]);
u2 = ValGradForward(u(2), [0 1 0]);
u3 = ValGradForward(u(3), [0 0 1]);


W = 0.5 * ( u1^2 + 2*u2^2 + u3^2 - 2*u1*u2 - 2*u2*u3) - (G(1) * u1 + G(2) * u2 + G(3) * u3);

AD = W.double;
val = AD(:,1);
grad = AD(:,2:end);
%grad(1,1) = 0;