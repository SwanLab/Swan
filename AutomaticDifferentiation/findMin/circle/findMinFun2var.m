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

figure(1)

plot(plotX); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotX) max(plotX)])

figure(2)

plot(plotY); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotY) max(plotY)])

figure(3)

plot(plotValX); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValX) max(plotValX)])

figure(4)

plot(plotValY); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValY) max(plotValY)])

valMin = x; % Value that minimizes the grad.
gradMin(1) = grad(1); % Gradient respect to x.
gradMin(2) = grad(2); % Gradient respect to y.
