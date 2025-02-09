function plottingMomentumParameter
n = 20000;
lambda(1) = 0;
gamma(1) = 0;
for i = 2:n
    lambda(i) = 0.5*(1 + sqrt(1+4*lambda(i-1)^2));
    gamma(i) = (lambda(i-1) - 1)/lambda(i);
end


plot(1:n,[lambda;gamma],'+-')
legend('\lambda','\gamma','Location','Best')
end