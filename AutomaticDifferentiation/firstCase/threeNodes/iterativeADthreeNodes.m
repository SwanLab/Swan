function [val, grad] = iterativeADthreeNodes(u0)

[~,numElem] = size(u0);

G = zeros(1,numElem);
G(1,numElem) = 1 * 10^(0);

A = 1 * 10^(-6);
E = 1 * 10^(9); %young modulus steel
L = 1 * 10^(3);

k = ((A*E)/L);

u(1) = ValGradForward(u0(1), [1 0 0]);
u(2) = ValGradForward(u0(2), [0 1 0]);
u(3) = ValGradForward(u0(3), [0 0 1]);

disp = (u(1)-u(2))^2;

for i = 2:(numElem-1)
    disp = disp + (-u(i-1)+2*u(i)-u(i+1))^2;
end

disp = disp + (-u(numElem-1)+u(numElem))^2;

disp = k^2 * disp;

force = G(1) * u(1);

for i = 2:numElem
    force = force + G(i) * u(i);
end

W = 0.5 * disp - force;

AD = W.double;
val = AD(:,1);
grad = AD(:,2:end);