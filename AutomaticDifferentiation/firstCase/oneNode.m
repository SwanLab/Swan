clear;
%% DEFINE INITIAL VARIABLES AND MATRIX
u0 = [1 1 1];
u = u0;
un = u0;

plotu = zeros(1000,3);
plotValu = zeros(1000,3);
gradMin = zeros(3,1);
iterations = 0;

%% ITERATE STARTING WITH A FOR WITH THE 3 VARIABLES

[val, grad] = iterativeADfirstCase(u0); %Initial Iteration

alpha = 0.1;

while abs(grad(1)) > 10^(-12) && abs(grad(2)) > 10^(-12) && abs(grad(3)) > 10^(-12) && iterations < 10^5 %while grad == 0 or iterations above 50
%while iterations < 10^4 %while grad == 0 or iterations above 50
    iterations = iterations + 1; %iterations counter

    [val, grad] = iterativeADfirstCase(u);

    un = u - alpha * grad;

    [valN, gradN] = iterativeADfirstCase(un);

    if iterations<100
        disp(val);
    end

    % if valN > val
    % 
    %     alpha = alpha / 2;
    % 
    % end

    u = un;

    plotu(iterations,1) = grad(1);
    plotu(iterations,2) = grad(2);
    plotu(iterations,3) = grad(3);

    plotValu(iterations,1) = u(1);
    plotValu(iterations,2) = u(2);
    plotValu(iterations,3) = u(3);

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



valMin = u; % Value that minimizes the grad.
gradMin(1) = grad(1); % Value that minimizes the grad.
gradMin(2) = grad(2); % Value that minimizes the grad.
gradMin(3) = grad(3); % Value that minimizes the grad.
