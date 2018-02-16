function Perimeter_analysis
% Reg
xmin = -0.5;
xmax = 0.5;


 %eta = epsilon/h

mu = 1;

n_eps = 6;
epsilon = 0.5.^[0:n_eps]';
he = ones(n_eps+1,1)*6.25*1e-4;
%he = [1*1e-1 5*1e-2 2.5*1e-2 1.25*1e-2 6.25*1e-3];
he2 = [2.5*1e-2  2.5*1e-2  2.5*1e-2 1.25*1e-2 6.25*1e-3 3.125*1e-3 1.5625*1e-3]';

%he = 0.5.^[8:n_eps+8];

%figure(3)
%plot(epsilon,'-+')
%hold on
%plot(1:length(epsilon),h*ones(length(epsilon),1),'r')
etas(:,1) = ones(n_eps+1,1);
etas(:,2) = 5*ones(n_eps+1,1);
etas(:,3) = 10*ones(n_eps+1,1);
etas(:,4) = epsilon./he;
etas(:,5) = epsilon./he2;


for ieta = 1:size(etas,2)

for iepsilon = 1:n_eps+1
    eps = epsilon(iepsilon);
   % h = he(iepsilon) 
    h = eps/etas(iepsilon,ieta);
    h_t(iepsilon,1) = h;
    N = round(1 + (xmax-xmin)/h);
    h = (xmax-xmin)/(N+1);
    x = linspace(xmin,xmax,N)';
    txi = compute_caracteristic(x);
    M = massMatrix(N);
    K = stiffMatrix(N);
    
    % Compute ve
    ve = (((mu*eps/h)^2)*K + M)\(M*txi);
    Pe(iepsilon) = (ve)'*(h*M)*(ve - 2*txi);
    vol = compute_volum(txi,h*M);
    Per(iepsilon,ieta) = 4/eps*[Pe(iepsilon) + vol]/mu;
    err(iepsilon,1) = ((ve - txi)'*(h*M)*(ve - txi))/((txi')*(h*M)*txi);

    figure(1)
    plot(Per,'-+')
    figure(2)
    plot(x,[txi,ve],'-+')
    figure(3)
    semilogy(err,'-+')
    figure(4)
    semilogy(epsilon(1:iepsilon),'-+')
%     figure(5)
%     loglog(h_t(1:iepsilon),err,'-+')    
%     eigK = eig(h*K);
%     Keig(:,iepsilon) = eigK(end-6:end);
%     eigM = eig(1/h*M);
%     Meig(:,iepsilon) = eigM(end-6:end);
    
    %condK = condest(h*K(2:end-1,2:end-1))
    %condM = condest(1/h*M(2:end-1,2:end-1))
    %axis([-0.5 0.5 0 1])

end


end


plot(epsilon,Per,'+')

xdata = log10(h_t);
ydata = log10(err);
f = fittype('a*x + b');
[fit1,gof,fitinfo] = fit(xdata,ydata,f,'StartPoint',[1 1]);
rate_of_conv = fit1.a;
Const = 10^fit1.b;

xdata = log10(err(end-5:end-1));
ydata = log10(err(end-4:end));
f = fittype('a*x + b');
[fit1,gof,fitinfo] = fit(xdata,ydata,f,'StartPoint',[1 1]);
rate_of_conv_num = fit1.a;
Const_num = 10^fit1.b;






end


function M = massMatrix(N)
D = sparse(1:N,1:N,4*ones(1,N),N,N);
E = sparse(2:N,1:N-1,ones(1,N-1),N,N);
M = E+D+E';
%M = diag(2*ones(N,1),0) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
M(1,1) = 2;
M(end,end) = 2;
M = 1/6*M;
end

function K = stiffMatrix(N)
D = sparse(1:N,1:N,2*ones(1,N),N,N);
E = sparse(2:N,1:N-1,-ones(1,N-1),N,N);
K = E+D+E';
K(1,1) = 1;
K(end,end) = 1;
end

function txi = compute_caracteristic(x)
txi = zeros(size(x));
txi(x > 0) = 1;
txi(x < 0) = 0;

%txi(x< -0.25 & x > 0.25) = 0;
%txi(x> -0.25 & x < 0.25) = 1;


end

function vol = compute_volum(txi,M)
vol = txi'*M*txi; 
end