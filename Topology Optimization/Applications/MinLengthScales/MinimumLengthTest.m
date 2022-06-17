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

% -----
% Choice A

rhoe = zeros(1,N);

for i=1:N % Loop over elements
    xcurr = 0.5*(x(i)+x(i+1));
    iN = find(abs(x-xcurr)<=R); % Neighbourhood
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

% Choice B
s.quadratureOrder = 'LINEAR';
s.mesh.nnodeElem  = 2;
s.mesh.nnodes = N+1;
s.mesh.type = 'LINE';
s.mesh.kFace = 0;
s.mesh.geometryType = 'Line';
s.mesh.coord = x';
for i=1:N
    connec(i,1) = i;
    connec(i,2) = i+1;
end
s.mesh.connec = connec;
s.mesh.nelem = N;
s.mesh.ndim = 1;
for i=1:N
    coordElem(:,:,i) = [x(i),x(i+1)];
end
s.mesh.coordElem = coordElem;
s.mesh.interpolation = Interpolation.create(s.mesh,'LINEAR');
Filter = Filter_P1_Density(s);
% -----

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