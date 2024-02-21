clear;

%% DEFINE INITIAL VARIABLES AND MATRIX

x0 = [1 1 1];
x = x0;
xn = x0;

plotX1 = zeros(1000,1);
plotX2 = zeros(1000,1);
plotX3 = zeros(1000,1);
plotValX1 = zeros(1000,1);
plotValX2 = zeros(1000,1);
plotValX3 = zeros(1000,1);
gradMin = zeros(3,1);
iterations = 0;

%% ITERATE STARTING WITH A FOR WITH THE 3 VARIABLES

[val, grad] = iterativeAD(x0); %Initial Iteration

alpha = 0.1;

while abs(grad(1)) > 10^(-12) && abs(grad(2)) > 10^(-12) && abs(grad(3)) > 10^(-12) || iterations < 10^4 %while grad == 0 or iterations above 50

    iterations = iterations + 1; %iterations counter

    [val, grad] = iterativeAD(x);

    xn = x - alpha * grad;

    [valN, gradN] = iterativeAD(xn);

    if valN > val

        alpha = alpha / 2;

    end

    x = xn;

    plotX1(iterations) = grad(1);
    plotX2(iterations) = grad(2);
    plotX3(iterations) = grad(3);

    plotValX1(iterations) = x(1);
    plotValX2(iterations) = x(2);
    plotValX3(iterations) = x(3);

end

%% PLOT

figure(1)

plot(plotX1); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotX1) max(plotX1)])

figure(2)

plot(plotX2); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotX2) max(plotX2)])

figure(3)

plot(plotX3); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotX3) max(plotX3)])

figure(4)

plot(plotValX1); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValX1) max(plotValX1)])

figure(5)

plot(plotValX2); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValX2) max(plotValX2)])

figure(6)

plot(plotValX3); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValX3) max(plotValX3)])



valMin = x; % Value that minimizes the grad.
gradMin(1) = grad(1); % Value that minimizes the grad.
gradMin(2) = grad(2); % Value that minimizes the grad.
gradMin(3) = grad(3); % Value that minimizes the grad.
