clear;
%% DEFINE INITIAL VARIABLES AND MATRIX
alpha = 0.001;
u = [1 1 1];

plotu = zeros(1000,3);
plotValu = zeros(1000,3);
plotVal = zeros(1000,1);
iterations = 0;

%% ITERATE STARTING WITH A FOR WITH THE 3 VARIABLES

[val(1), grad(1,:)] = iterativeADOneNodeNewhookean(u); %Initial Iteration

u = u - alpha * grad(1,:);

while abs(grad(1)) > 10^(-12) && abs(grad(2)) > 10^(-12) && abs(grad(3)) > 10^(-12) && iterations < 10^4 %while grad == 0 or iterations above 50
    iterations = iterations + 1; %iterations counter

    [val(2), grad(2,:)] = iterativeADOneNodeNewhookean(u);

    u = u - alpha * grad(2,:);

    val(1) = val(2);
    grad(1,:) = grad(2,:);

    % if val(2) > val(1)
    %
    %     alpha = alpha / 2;
    %
    % end

    plotu(iterations,1) = grad(1,1);
    plotu(iterations,2) = grad(1,2);
    plotu(iterations,3) = grad(1,3);

    plotValu(iterations,1) = u(1);
    plotValu(iterations,2) = u(2);
    plotValu(iterations,3) = u(3);

    plotVal(iterations) = val(1);

end

%% PLOT

figure(1)

plot(plotu(:,1)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotu(:,1)) max(plotu(:,1))])

figure(2)

plot(plotu(:,2)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotu(:,2)) max(plotu(:,2))])

figure(3)

plot(plotu(:,3)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotu(:,3)) max(plotu(:,3))])

figure(4)

plot(plotValu(:,1)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValu(:,1)) max(plotValu(:,1))])

figure(5)

plot(plotValu(:,2)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValu(:,2)) max(plotValu(:,2))])

figure(6)

plot(plotValu(:,3)); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotValu(:,3)) max(plotValu(:,3))])

figure(7)

plot(plotVal); %plot of the gradient tending to the min.

xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotVal) max(plotVal)])


valMin = u; % Value that minimizes the grad.
gradMin(1) = grad(1,1); % Value that minimizes the grad.
gradMin(2) = grad(1,2); % Value that minimizes the grad.
gradMin(3) = grad(1,3); % Value that minimizes the grad.