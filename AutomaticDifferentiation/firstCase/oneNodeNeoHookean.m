clear;
%% DEFINE INITIAL VARIABLES AND MATRIX
alpha = 0.001;
u = [1 1 1];

plotGrad = zeros(100,3);
plotU = zeros(100,3);
plotVal = zeros(100,1);
iterations = 0;

%% ITERATE STARTING WITH A FOR WITH THE 3 VARIABLES

[val(1), grad(1,:)] = iterativeADOneNodeNeoHookean(u); %Initial Iteration

u = u - alpha * grad(1,:);

while abs(grad(1)) > 10^(-12) && abs(grad(2)) > 10^(-12) && abs(grad(3)) > 10^(-12) && iterations < 10^4 %while grad == 0 or iterations above 50
    iterations = iterations + 1; %iterations counter

    [val(2), grad(2,:)] = iterativeADOneNodeNeoHookean(u);

    u = u - alpha * grad(2,:);

    val(1) = val(2);
    grad(1,:) = grad(2,:);

    if val(2) > val(1)

        alpha = alpha / 2;

    end

    plotGrad(iterations,1) = grad(1,1);
    plotGrad(iterations,2) = grad(1,2);
    plotGrad(iterations,3) = grad(1,3);

    plotU(iterations,1) = u(1);
    plotU(iterations,2) = u(2);
    plotU(iterations,3) = u(3);

    plotVal(iterations) = val(1);

end

%% PLOT

figure(1)

plot(plotGrad(:,1)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotGrad(:,1)) max(plotGrad(:,1))])

figure(2)

plot(plotGrad(:,2)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotGrad(:,2)) max(plotGrad(:,2))])

figure(3)

plot(plotGrad(:,3)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotGrad(:,3)) max(plotGrad(:,3))])

figure(4)

plot(plotU(:,1)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,1)) max(plotU(:,1))])

figure(5)

plot(plotU(:,2)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,2)) max(plotU(:,2))])

figure(6)

plot(plotU(:,3)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,3)) max(plotU(:,3))])

figure(7)

plot(plotVal); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotVal) max(plotVal)])


valMin = u; % Value that minimizes the grad.
gradMin(1) = grad(1,1); % Value that minimizes the grad.
gradMin(2) = grad(1,2); % Value that minimizes the grad.
gradMin(3) = grad(1,3); % Value that minimizes the grad.