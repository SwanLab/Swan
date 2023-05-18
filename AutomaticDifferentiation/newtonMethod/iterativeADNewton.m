function [val, grad, grad2] = iterativeADNewton(u)

G = [0; 0; 1];
G = transpose(G);

u1 = ValGradForward(ValGradForward(0, [1 0 0]),1);
u2 = ValGradForward(ValGradForward(u(2), [0 1 0]),1);
u3 = ValGradForward(ValGradForward(u(3), [0 0 1]),1);


W = 0.5 * ( u1^2 + 2*u2^2 + u3^2 - 2*u1*u2 - 2*u2*u3) - (G(1) * u1 + G(2) * u2 + G(3) * u3);

firstGrad = W.val.double;
secondGrad = W.grad.double;

val = firstGrad(:,1);
grad = firstGrad(:,2:end);
grad2 = secondGrad(:,2:end);