%% Minimum length scales test

clear
clc

N   = 500; % Number of elements
x   = linspace(0,1,N+1); % Number of nodes

% Density - initial map example
rho = zeros(1,round(0.25*length(x)));
rho = [rho,ones(1,round(0.125*length(x)))];
rho = [rho,zeros(1,round(0.25*length(x)))];
rho = [rho,ones(1,round(0.125*length(x)))];
rho = [rho,zeros(1,length(x)-length(rho))];

R = 0.10*max(x); % Filter radius

for i=1:N % Loop over elements
    xcurr = 0.5*(x(i)+x(i+1));
    iN = find(abs(x-xcurr)<=R); % Neighborhood
    Num = 0;
    Den = 0;
    for j=1:length(iN)-1
        jN = iN(j);
        Vj = x(jN+1)-x(jN);
        xj = 0.5*(x(jN)+x(jN+1));
        rhoj = 0.5*(rho(jN+1)+rho(jN));
        Num = Num + (R-abs(xcurr-xj))*Vj*rhoj;
        Den = Den + (R-abs(xcurr-xj))*Vj;
    end
    rhoe(i) = Num/Den;
end

xc = 0.5*(x(1:end-1)+x(2:end));
rhoe = interp1(xc,rhoe,x);

beta=120;
eta=0;
rhoP=(tanh(beta*eta)+tanh(beta*(rhoe-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta)));

figure
plot(x,rho,'-k',x,rhoe,x,rhoP,'--k')
grid on
grid minor
xlabel('x coor')
ylabel('density')
legend('Real','Regularized','Projected')
ylim([-0.1 1.1])

