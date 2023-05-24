function [plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinNVariables(u)

iterations = 1;

plotU = zeros(iterations,3);
plotVal = zeros(iterations,1);
plotGrad = zeros(iterations,1);
plotGrad2 = zeros(iterations,1);

iterations = 0;

grad = 1;

while abs(grad) > 10^(-6)
    iterations = iterations + 1;
    [val, grad, grad2] = iterativeADNewton(u);
    h = grad2.^(-1);
    h(isinf(h)) = 0;
    h = [h(1,1) h(2,2) h(3,3)];
    % u = u - transpose(h * transpose(grad));
    u = u - h .* grad;

    plotU(iterations,1:3) = u;
    plotVal(iterations) = val;
    plotGrad(iterations,1:3) = grad;
    plotGrad2(iterations,1:3) = [grad2(1,1),grad2(2,2),grad2(3,3)];
end

%figure(1); plot(plotU(:,1)); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,1)) max(plotU(:,1))]);
figure(2); plot(plotU(:,2)); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,2)) max(plotU(:,2))]);
figure(3); plot(plotU(:,3)); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,3)) max(plotU(:,3))]);
%figure(4); plot(plotGrad(:,1)); xlabel("Num. of iterations"); ylabel("First Gradient"); grid; axis([0 iterations min(plotGrad(:,1)) max(plotGrad(:,1))]);
figure(5); plot(plotGrad(:,2)); xlabel("Num. of iterations"); ylabel("First Gradient"); grid; axis([0 iterations min(plotGrad(:,2)) max(plotGrad(:,2))]);
figure(6); plot(plotGrad(:,3)); xlabel("Num. of iterations"); ylabel("First Gradient"); grid; axis([0 iterations min(plotGrad(:,3)) max(plotGrad(:,3))]);
%figure(7); plot(plotGrad2(:,1)); xlabel("Num. of iterations"); ylabel("Second Gradient"); grid; axis([0 iterations min(plotGrad2(:,1)) max(plotGrad2(:,1))]);
%figure(8); plot(plotGrad2(:,2)); xlabel("Num. of iterations"); ylabel("Second Gradient"); grid; axis([0 iterations min(plotGrad2(:,2)) max(plotGrad2(:,2))]);
%figure(9); plot(plotGrad2(:,3)); xlabel("Num. of iterations"); ylabel("Second Gradient"); grid; axis([0 iterations min(plotGrad2(:,3)) max(plotGrad2(:,3))]);
figure(10); plot(plotVal); xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotVal) max(plotVal)]);