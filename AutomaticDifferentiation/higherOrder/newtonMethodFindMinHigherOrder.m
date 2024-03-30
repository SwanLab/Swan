function [plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinHigherOrder(u)
iterations = 11854; %41347
plotU = zeros(iterations,1);
plotVal = zeros(iterations,1);
plotGrad = zeros(iterations,1);
plotGrad2 = zeros(iterations,1);
alpha = 0.0001;
iterations = 0;
grad = 1;
while abs(grad) > 10^(-6)
    iterations = iterations + 1;
    [val, grad, grad2] = iterativeADHigherOrder(u);
    u = u - grad2^(-1) * grad;
    plotU(iterations) = u;
    plotVal(iterations) = val;
    plotGrad(iterations) = grad;
    plotGrad2(iterations) = grad2;
end

figure(1);
plot(plotU); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU) max(plotU)]);
figure(2);
plot(plotVal); xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotVal) max(plotVal)]);
figure(3);
plot(plotGrad); xlabel("Num. of iterations"); ylabel("First Gradient"); grid; axis([0 iterations min(plotGrad) max(plotGrad)]);
figure(4);
plot(plotGrad2); xlabel("Num. of iterations"); ylabel("Second Gradient"); grid; axis([0 iterations min(plotGrad2) max(plotGrad2)]);
figure(5);
plot(plotU,plotVal); xlabel("X values"); ylabel("Function Value"); title("Part of the function evaluated through iterations."); grid; axis([min(plotU) max(plotU) min(plotVal) max(plotVal)]);