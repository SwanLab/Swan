clear;
plotf = zeros(100,1);
plotff = zeros(100,1);
plotfff = zeros(100,1);
y = -50:1:49;

for i = 1 : 100

    x = ValGradForward(ValGradForward(-51+i,1),1);
    f = x^3;
    plotf(i,1) = f.val.val;
    plotff(i,1) = f.grad.val;
    plotfff(i,1) = f.grad.grad;

end

figure(1);
plot(y,plotf);
figure(2);
plot(y,plotff);
figure(3);
plot(y,plotfff);