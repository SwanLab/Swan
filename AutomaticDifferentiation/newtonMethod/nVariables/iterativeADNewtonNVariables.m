function [val, grad, grad2] = iterativeADNewtonNVariables(u0)

[~,numElem] = size(u0);

%u = zeros(1,numElem);
%disp = 0;

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


%W = 0.5 * ( u(1)^2 + 2*u(2)^2 + u(3)^2 - 2*u(1)*u(2) - 2*u(2)*u(3)) - (G(1) * u(1) + G(2) * u(2) + G(3) * u(3));

disp = (u(1)-u(2))^2;

for i = 2:(numElem-1)
    disp = disp + (u(i)-u(i+1))^2;
end


force = G(1) * u(1);

for i = 2:numElem
    force = force + G(i) * u(i);
end

W = 0.5 * disp - force;

firstGrad = W.val.double;
secondGrad = W.grad.double;

val = firstGrad(:,1);
grad = firstGrad(:,2:end);
grad2 = secondGrad(1:end,(1+numElem):end);