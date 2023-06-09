clear;

%% DEFINE INITIAL VARIABLES AND MATRIX

x = [0 0];

plotX = zeros(50,1);
plotY = zeros(50,1);
plotValX = zeros(50,1);
plotValY = zeros(50,1);

iterations = 0;

%% ITERATE STARTING WITH A FOR WITH THE 2 VARIABLES

[~,numElem] = size(x);

grad = 1 + zeros(1,numElem);

alpha = 0.1;

while abs(sum(grad(1,1:end))) > 10^(-10) && iterations < 10^5

    iterations = iterations + 1; %iterations counter

    [val, grad] = iterativeAD2var(x);

    x = x - alpha * grad;

    %% PLOT

    plotX(iterations) = grad(1);
    plotY(iterations) = grad(2);

    plotValX(iterations) = x(1);
    plotValY(iterations) = x(2);


end

% Create plots
t = tiledlayout(2,2);

nexttile

plot(plotX); xlabel("Num. of iterations"); ylabel("Gradient"); title("Gradient on direction X"); grid; axis([0 iterations -1 1]) %plot of the gradient tending to the min. 

nexttile

plot(plotY); xlabel("Num. of iterations"); ylabel("Gradient"); title("Gradient on direction Y"); grid; axis([0 iterations -1 1]) %plot of the gradient tending to the min.

nexttile

plot(plotValX); xlabel("Num. of iterations"); ylabel("Value"); title("X values"); grid; axis([0 iterations 0 2]) %plot of the gradient tending to the min.

nexttile

plot(plotValY); xlabel("Num. of iterations"); ylabel("Value"); title("Y values"); grid; axis([0 iterations 0 2]) %plot of the gradient tending to the min.



valMin = x; % Value that minimizes the grad.
gradMin(1) = grad(1); % Gradient respect to x.
gradMin(2) = grad(2); % Gradient respect to y.
