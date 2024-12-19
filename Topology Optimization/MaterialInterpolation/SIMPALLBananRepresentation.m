function SIMPALL(E1,E0,nu1,nu0, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMPALL(1,0.01,1/3,1/3, N)         A %%      +3
%% SIMPALL(1,0.001, 0, 1/3, N)        B %%      +3
%% SIMPALL(1,0.01, -0.5, 1/3, N)      C %%       4
%% SIMPALL(1,0.001, -0.75, 1/3, N)    D %%       3
%% SIMPALL(1,0.001, -0.9, 1/3, N)     E %%       3
%% SIMPALL(1,1/3, 0.3, 0.35, N)      bi %%      +3
%% SIMPALL(1,0.9, 0.3, -0.9, N)      bi %%       3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nota la dimensio maxima imposable depen de Poisson (nu) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=@(E,nu) E/(2*(1+nu));
kappa = @(E,nu, N) E/(N*(1+nu-N*nu));
mu0 = mu(E0,nu0);
mu1 = mu(E1,nu1); 
kappa0 = kappa(E0, nu0, N);
kappa1 = kappa(E1, nu1, N);

etaMu = @(mu, kappa, N)  -(mu*(4*mu - kappa*N^2 - 2*mu*N^2 + 2*mu*N))/(2*N*(kappa + 2*mu)); 
etaKappa = @(mu, kappa, N) (2*mu*(N - 1))/N; 

% SIMP-ALL coeficients
A = @(f0, f1, df0, df1) -(-f0^2+2*f0*f1-f1^2+df0*df1) /(df1-f1+f0);
B = @(f0, f1, df0, df1) -(df0*f1-df0*df1+df1*f0-2*f0*f1+2*f0^2)/(df1-f1+f0);
C = @(f0)                f0;
D = @(f0, f1, df0, df1) -(df0+df1+2*f0-2*f1)/(df1-f1+f0);


df0 = @(eta, f0, f1) (f0+eta)*(f1-f0)/(f1+eta);
df1 = @(eta, f0, f1) (f1+eta)*(f1-f0)/(f0+eta);


%% Case mu
F0 = mu0;
F1 = mu1;
Eta0 = etaMu(mu0, kappa0, N);
Eta1 = etaMu(mu1, kappa1, N);
DF0 = df0(Eta0, F0, F1);
DF1 = df1(Eta1, F0, F1);
data = computeData(F0, F1, DF0, DF1, Eta0, Eta1, A, B, C, D);
computePlot('Shear modulus', data);

%% Case kappa
F0 = kappa0;
F1 = kappa1;
Eta0 = etaKappa(mu0, kappa0, N);
Eta1 = etaKappa(mu1, kappa1, N);
DF0 = df0(Eta0, F0, F1);
DF1 = df1(Eta1, F0, F1);
data = computeData(F0, F1, DF0, DF1, Eta0, Eta1, A, B, C, D);
computePlot('Bulk modulus', data);

end

function computePlot(Title, data)
AA = data.AA;
BB = data.BB;
CC = data.CC;
DD = data.DD;
F0 = data.F0;
F1 = data.F1;
Eta0 = data.Eta0;
Eta1 = data.Eta1;
figure;
hold on
f = @(rho, A, B, C, D) (A.*rho.^2+B.*rho+C)./(D.*rho+1); %SIMPALL
plot(linspace(0,1,100),f(linspace(0,1,100), AA, BB, CC, DD), 'r', 'LineWidth',2); %SIMP-ALL plot
SIMP = @(rho) (1-rho.^3).*F0+rho.^3.*F1; %SIMP
plot(linspace(0,1,100),SIMP(linspace(0,1,100)), 'b', 'LineWidth',2); %SIMP plot
fHs = @(rho, Eta) F0.*(1-rho)+F1.*rho-((1-rho).*rho*(F1-F0)^2)./(F0.*rho+F1.*(1-rho)+Eta);
plot(linspace(0,1,100),fHs(linspace(0,1,100), Eta0), 'k', 'LineWidth',2); %LB plot
plot(linspace(0,1,100),fHs(linspace(0,1,100), Eta1), 'k', 'LineWidth',2); %UB plot
set(gca,'fontsize',16)
set(gca,'XTick',0:0.5:1)
title(Title,'fontsize',20)
xlabel('Density','fontsize',20)
legend({'SIMP-ALL','SIMP', 'HS-LB','HS-UB'},'location','best','fontsize',20)
daspect([1 1 2])
end

function data = computeData(F0, F1, DF0, DF1, Eta0, Eta1, A, B, C, D)
data.AA = A(F0, F1, DF0, DF1);
data.BB = B(F0, F1, DF0, DF1);
data.CC = C(F0);
data.DD = D(F0, F1, DF0, DF1);
data.F0 = F0;
data.F1 = F1;
data.Eta0 = Eta0;
data.Eta1 = Eta1;
end