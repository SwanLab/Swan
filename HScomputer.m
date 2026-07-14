clc,clear,close all
%% General
E  = 1;
nuVal = 1/3;
k2D  = @(nu) E./(2.*(1-nu));
mu2D = @(nu) E./(2.*(1+nu));
k3D  = @(nu) E./(3 -6.*nu);
mu3D = @(nu) E./(2.*(1+nu));
l3D  = @(nu) E*nu/((1+nu)*(1-2*nu));

%% H-S 2D
% Preliminaries
k0     = 1e-10; 
mu0    = 1e-10;
k1     = @(nu) k2D(nu);
mu1    = @(nu) mu2D(nu);

etak02D  = mu0;
etamu02D = (k0.*mu0)./(2.*mu0+k0);
etak12D  = @(nu) mu1(nu);
etamu12D = @(nu) (k1(nu).*mu1(nu))./(2.*mu1(nu)+k1(nu));

% H-S definition
f02D = k0; 
f12D = @(nu) k1(nu);
eta2D = @(nu) etak12D(nu);

num = @(rho,nu) (1-rho).*rho.*(f12D(nu) - f02D).^2;
den =  @(rho,nu) f02D.*rho + f12D(nu).*(1-rho) + eta2D(nu);
fHS2D = @(rho,nu) f02D.*(1-rho) + f12D(nu).*rho - num(rho,nu)./den(rho,nu); 

fHSdmg = @(phi,nu) fHS2D(1-phi,nu);
kUB2D  = @(phi,nu) (k0.*(phi) + k1(nu).*(1-phi) - ((1-phi).*phi.*(k1(nu)-k0).^2)./(k0.*(1-phi) + k1(nu).*phi + etak12D(nu)));

figure(1)
hold on
fplot(@(rho) fHS2D(rho,nuVal),[0 1])
fplot(@(phi) fHSdmg(phi,nuVal),[0 1])
fplot(@(phi) kUB2D(phi,nuVal),[0 1],'LineStyle','--')
legend('H-S 2D (density)','H-S 2D (damage)','kUB 2D')

%% H-S 2D bulk nD
% Preliminaries
k1n = @(nu,N) (N-2).*k3D(nu) - (N-3).*k2D(nu);
mu1n = @(nu,N) (N-2).*mu3D(nu) - (N-3).*mu2D(nu); % Not necessary. same 2D/3D

% H-S definition
den1 = @(nu,N) 1./(k0-k1n(nu,N));
den2 = @(rho,nu,N) rho./(k1n(nu,N) + 2.*mu1n(nu,N) - 2.*mu1n(nu,N)./N);
kUBnD = @(rho,nu,N) k1n(nu,N) + (1-rho)./(den1(nu,N)+den2(rho,nu,N));

kUBnDdmg = @(phi,nu,N) kUBnD(1-phi,nu,N);

figure(2)
hold on
fplot(@(rho) fHS2D(rho,nuVal),[0 1])
fplot(@(rho) kUBnD(rho,nuVal,2),[0 1],'LineStyle','--')
fplot(@(phi) kUBnDdmg(phi,nuVal,2),[0 1])
fplot(@(phi) kUB2D(phi,nuVal),[0 1],'LineStyle','--')
legend('H-S 2D (density)','kUB nD/2D (density)','kUB nD/2D (damage)','kUB 2D')

%% H-S ND (3D)
% Preliminaries
etak03D  = @(N) 2.*mu0.*(N-1)/N;
etamu03D = @(N) (mu0./(2.*mu0 + k0)).*( (mu0.*((N-2).*(N-1))./N) + k0.*N/2 );
etak13D  = @(nu,N) 2.*mu1n(nu,N).*(N-1)/N;
etamu13D = @(nu,N) (mu1n(nu,N)./(2.*mu1n(nu,N) + k1n(nu,N))).*( (mu1n(nu,N).*((N-2).*(N+1))./N) + k1n(nu,N).*N/2 );

% H-S definition
f03D = k0;
f13D = @(nu,N) k1n(nu,N);
eta3D = @(nu,N) etak13D(nu,N);

num3D = @(rho,nu,N) (1-rho).*rho.*(f13D(nu,N) - f03D).^2;
den3D = @(rho,nu,N) eta3D(nu,N)./(f03D.*f13D(nu,N)) + (rho./f13D(nu,N) + (1-rho)./f03D);
fHS3D = @(rho,nu,N) f03D.*(1-rho) + f13D(nu,N).*rho - (1./(f03D.*f13D(nu,N))).*(num3D(rho,nu,N)./den3D(rho,nu,N));

fHS3Ddmg = @(phi,nu,N) fHS3D(1-phi,nu,N);

figure(3)
hold on
fplot(@(rho) fHS3D(rho,nuVal,3),[0 1])
fplot(@(rho) kUBnD(rho,nuVal,3),[0 1],'LineStyle','--')
fplot(@(phi) fHS3Ddmg(phi,nuVal,3),[0 1])
fplot(@(phi) kUBnDdmg(phi,nuVal,3),[0 1],'LineStyle','--')
legend('H-S nD/3D (density)','kUB nD/3D (density)','H-S nD/3D (damage)','kUB nD/3D (damage)')

%% H-S ND (2D)
figure(4)
hold on
fplot(@(rho) fHS3D(rho,nuVal,2),[0 1])
fplot(@(rho) kUBnD(rho,nuVal,2),[0 1],'LineStyle','--')
fplot(@(phi) fHS3Ddmg(phi,nuVal,2),[0 1])
fplot(@(phi) kUBnDdmg(phi,nuVal,2),[0 1],'LineStyle','--')
legend('H-S nD/2D (density)','kUB nD/2D (density)','H-S nD/2D (damage)','kUB nD/2D (damage)')


%% Plots verification
kUB3D   = @(phi,nu) k3D(nu).*(1-phi).*((2.*mu3D(nu)+l3D(nu)-k3D(nu))./(2.*mu3D(nu)+l3D(nu)-k3D(nu).*(1-phi)));
figure(5)
hold on
fplot(@(phi) fHS3Ddmg(phi,nuVal,3),[0 1])
fplot(@(phi) kUBnDdmg(phi,nuVal,3),[0 1],'LineStyle','--')
fplot(@(phi) kUB3D(phi,nuVal),[0 1],'LineStyle','--')
legend('H-S nD/3D','kUB nD/3D','kUB 3D')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f02D = k0; 
f12D = @(nu) mu1(nu);
eta2D = @(nu) etamu12D(nu);

num = @(rho,nu) (1-rho).*rho.*(f12D(nu) - f02D).^2;
den =  @(rho,nu) f02D.*rho + f12D(nu).*(1-rho) + eta2D(nu);
fHS2D = @(rho,nu) f02D.*(1-rho) + f12D(nu).*rho - num(rho,nu)./den(rho,nu); 

fHSdmg = @(phi,nu) fHS2D(1-phi,nu);
muUB2D = @(phi,nu) (mu0.*(phi) + mu1(nu).*(1-phi) - ((1-phi).*phi.*(mu1(nu)-mu0).^2)./(mu0.*(1-phi) + mu1(nu).*phi + etamu12D(nu)));

figure(6)
hold on
fplot(@(rho) fHS2D(rho,nuVal),[0 1])
fplot(@(phi) fHSdmg(phi,nuVal),[0 1])
fplot(@(phi) muUB2D(phi,nuVal),[0 1],'LineStyle','--')
legend('H-S 2D (density)','H-S 2D (damage)','muUB 2D')

%%
den1 = @(nu,N) 1./(mu0-mu1n(nu,N));
dennum = @(rho,nu,N) 2.*N.*rho.*(k1n(nu,N) + 2.*mu1n(nu,N)).*(N-1);
denden = @(rho,nu,N) mu1n(nu,N).*(N^2 +N-2).*(k1n(nu,N).*N + 2.*mu1n(nu,N).*(N-1));
den2 = @(rho,nu,N) dennum(rho,nu,N)./denden(rho,nu,N);
muUBnD = @(rho,nu,N) mu1n(nu,N) + (1-rho)./(den1(nu,N)+den2(rho,nu,N));

muUBnDdmg = @(phi,nu,N) muUBnD(1-phi,nu,N);

figure(7)
hold on
fplot(@(rho) fHS2D(rho,nuVal),[0 1])
fplot(@(rho) muUBnD(rho,nuVal,2),[0 1],'LineStyle','--')
fplot(@(phi) muUBnDdmg(phi,nuVal,2),[0 1])
fplot(@(phi) muUB2D(phi,nuVal),[0 1],'LineStyle','--')
legend('H-S 2D (density)','muUB nD/2D (density)','muUB nD/2D (damage)','muUB 2D')

%%
f03D = mu0;
f13D = @(nu,N) mu1n(nu,N);
eta3D = @(nu,N) etamu13D(nu,N);

num3D = @(rho,nu,N) (1-rho).*rho.*(f13D(nu,N) - f03D).^2;
den3D = @(rho,nu,N) eta3D(nu,N)./(f03D.*f13D(nu,N)) + (rho./f13D(nu,N) + (1-rho)./f03D);
fHS3D = @(rho,nu,N) f03D.*(1-rho) + f13D(nu,N).*rho - (1./(f03D.*f13D(nu,N))).*(num3D(rho,nu,N)./den3D(rho,nu,N));

fHS3Ddmg = @(phi,nu,N) fHS3D(1-phi,nu,N);

figure(8)
hold on
fplot(@(rho) fHS3D(rho,nuVal,3),[0 1])
fplot(@(rho) muUBnD(rho,nuVal,3),[0 1],'LineStyle','--')
fplot(@(phi) fHS3Ddmg(phi,nuVal,3),[0 1])
fplot(@(phi) muUBnDdmg(phi,nuVal,3),[0 1],'LineStyle','--')
legend('H-S nD/3D (density)','muUB nD/3D (density)','H-S nD/3D (damage)','muUB nD/3D (damage)')

%% H-S ND (2D)
figure(9)
hold on
fplot(@(rho) fHS3D(rho,nuVal,2),[0 1])
fplot(@(rho) muUBnD(rho,nuVal,2),[0 1],'LineStyle','--')
fplot(@(phi) fHS3Ddmg(phi,nuVal,2),[0 1])
fplot(@(phi) muUBnDdmg(phi,nuVal,2),[0 1],'LineStyle','--')
legend('H-S nD/2D (density)','muUB nD/2D (density)','H-S nD/2D (damage)','muUB nD/2D (damage)')


%% Plots verification
muUB3D  = @(phi,nu) mu3D(nu).*(1-phi).*((10.*(2.*mu3D(nu)+l3D(nu)) - 4.*(k3D(nu)+2.*mu3D(nu)))./(10.*(2.*mu3D(nu)+l3D(nu))-(1-phi).*4.*(k3D(nu)+2.*mu3D(nu))));

figure(10)
hold on
fplot(@(phi) fHS3Ddmg(phi,nuVal,3),[0 1])
fplot(@(phi) muUBnDdmg(phi,nuVal,3),[0 1])
fplot(@(phi) muUB3D(phi,nuVal),[0 1],'LineStyle','--')
legend('H-S nD/3D','muUB nD/3D','muUB 3D')
