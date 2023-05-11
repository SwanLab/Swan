clear;

%% DEFINE INITIAL VARIABLES AND MATRIX

x0 = [1 1 1];
x = x0;
xn = x0;

plotX = zeros(1000,3);
valMin = zeros(3,1);
gradMin = zeros(3,1);
iterations = zeros(3,1);

%% ITERATE STARTING WITH A FOR WITH THE 3 VARIABLES

for i = 1:3

    [val, grad] = iterativeAD(x0); %Initial Iteration

    alpha = 1;

    while abs(grad(i)) > 10^(-12) && iterations(i) < 50 %while grad == 0 or iterations above 50

        iterations(i) = iterations(i) + 1; %iterations counter

        [val, grad] = iterativeAD(x);

        xn(i) = x(i) - alpha * grad(i);

        [valN, gradN] = iterativeAD(xn);

        if valN > val

            alpha = alpha / 2;

        end

        x(i) = xn(i);

        %% PLOT

        plotX(iterations(i),i) = grad(i);

        figure(i)

        plot(plotX(:,i)); %plot of the gradient tending to the min.

        xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 50 -1 1])


    end

    valMin(i) = val; % Value that minimizes the grad.
    gradMin(i) = grad(i); % Value that minimizes the grad.

end
