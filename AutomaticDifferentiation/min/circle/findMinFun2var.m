clear;

%% DEFINE INITIAL VARIABLES AND MATRIX

x0 = [0 0];
x = x0;
xn = x0;

plotX = zeros(50,1);
plotY = zeros(50,1);
plotValX = zeros(50,1);
plotValY = zeros(50,1);
gradMin = zeros(2,1);
iterations = 0;

%% ITERATE STARTING WITH A FOR WITH THE 2 VARIABLES

[val, grad] = iterativeAD2var(x0); %Initial Iteration

alpha = 1;

while abs(grad(1)) > 10^(-12) && abs(grad(2)) > 10^(-12) && iterations < 10^4 %while grad == 0 or iterations above 50

    iterations = iterations + 1; %iterations counter

    [val, grad] = iterativeAD2var(x);

    xn = x - alpha * grad;

    [valN, gradN] = iterativeAD2var(xn);

    if valN > val

        alpha = alpha / 2;

    end

    x = xn;

    %% PLOT

    plotX(iterations) = grad(1);
    plotY(iterations) = grad(2);

    plotValX(iterations) = x(1);
    plotValY(iterations) = x(2);


end

figure(1)

plot(plotX); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations -1 1])

figure(2)

plot(plotY); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations -1 1])

figure(3)

plot(plotValX); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations -1 1])

figure(4)

plot(plotValY); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations -1 1])

valMin = x; % Value that minimizes the grad.
gradMin(1) = grad(1); % Gradient respect to x.
gradMin(2) = grad(2); % Gradient respect to y.
