clear;
nodes = 50;
u0 = 1 + zeros(1,nodes);
[plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinNVariables(u0);


load('C:\Users\usuari\OneDrive\Escritorio\TFG\Swan\AutomaticDifferentiation\newtonMethod\nVariables\iterationsVSnodes.mat')
i = 3:50;
j = 1:27287;

gradquad = plotGrad(:,2:end).^2;
gradquad1 = sum(gradquad,2);
norma = sqrt(gradquad1);

figure(1); loglog(i,time); xlabel("Log of nodes"); ylabel("Log of iterations"); grid;
figure(2); semilogx(norma,j); xlabel("Log norma gradient"); ylabel("Iterations"); grid;

