clear;

%% DEFINE INITIAL VARIABLES AND MATRIX
u = [1 1 1];

iterations = 349;

plotGrad = zeros(iterations,3);
plotU = zeros(iterations,1);
plotVal = zeros(iterations,1);

iterations = 0;

%% ITERATE STARTING WITH A FOR WITH THE 2 VARIABLES

[~,numElem] = size(u);

grad = 1 + zeros(1,numElem);

alpha = 1 * 10^(-2);
tic
while abs(sum(grad(1,2:end))) > 10^(-6) && iterations < 10^3 %while grad == 0 or iterations above 50
    %while iterations < 10^4 %while grad == 0 or iterations above 50
    iterations = iterations + 1; %iterations counter
    [val, grad] = iterativeADThreeNodeNeoHookean_REVERSE(u);
    [val, grad] = iterativeADThreeNodeNeoHookean(u);

    u = u - alpha * grad;

    %% PLOT

    plotGrad(iterations,1) = grad(1,1);
    plotGrad(iterations,2) = grad(1,2);
    plotGrad(iterations,3) = grad(1,3);

    plotU(iterations,1) = u(1);
    plotU(iterations,2) = u(2);
    plotU(iterations,3) = u(3);

    plotVal(iterations) = val;


end
toc
%% PLOT

t = tiledlayout(2,2);

nexttile

plot(plotGrad(:,2)); xlabel("Num. of iterations"); ylabel("Gradient"); title("Gradient on node 2"); grid; axis([0 iterations -0.3 0.2]) %plot of the gradient tending to the min.

nexttile

plot(plotGrad(:,3)); xlabel("Num. of iterations"); ylabel("Gradient"); title("Gradient on node 3"); grid; axis([0 iterations -0.3 0.2]) %plot of the gradient tending to the min.

nexttile

plot(plotU(:,2)); xlabel("Num. of iterations"); ylabel("Displacement"); title("Displacement on node 2"); grid; axis([0 iterations 0 2]) %plot of the gradient tending to the min.

nexttile

plot(plotU(:,3)); xlabel("Num. of iterations"); ylabel("Displacement"); title("Displacement on node 3"); grid; axis([0 iterations 0 2]) %plot of the gradient tending to the min.