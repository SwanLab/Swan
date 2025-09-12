load('cost.mat')
load('constraint.mat')

figure;

yyaxis left
plot(iterations, -cost, '-.')
ylabel('First Eigenvalue')
yyaxis right
plot(iterations, 0.4*(constraint+1), '-.')
ylabel('Constraint')
xlabel('Iteration')
grid
xlabel('Iterations')