load('HomogenizationResults3.mat');   
% Parâmetros SIMP
E0   = 1.0;
nu   = 0.30;
p    = 3;
Emin = 1e-6;
Efun = @(phi) Emin + (E0 - Emin).*phi.^p;


lam_psn = @(phi) Efun(phi).*nu./((1+nu).*(1-2*nu));
mu_psn  = @(phi) Efun(phi)./(2*(1+nu));
C1111_simp_raw = @(phi) lam_psn(phi) + 2.*mu_psn(phi);


scale_C1111 = 1.1538 / C1111_simp_raw(1);
C1111_simp  = @(phi) scale_C1111 * C1111_simp_raw(phi);   


C1122_simp = @(phi) Efun(phi).*nu./(1 - nu^2);      


C1212_simp = @(phi) 2.*mu(phi);   



% --------- Plots ----------
tiledlayout(1,3)


nexttile; hold on; grid on
plot(volFrac, squeeze(Chomog(1,1,1,1,:)), 'o', 'DisplayName','Homog C_{1111}')
if exist('Interpolation','var'); fplot(Interpolation.fun(1,1,1,1), [0 1], 'DisplayName','Fit/Interp'); end
fplot(C1111_simp, [0 1], 'DisplayName','SIMP')
xlabel('\phi'); ylabel('C_{1111}'); legend


nexttile; hold on; grid on
plot(volFrac, squeeze(Chomog(1,1,2,2,:)), 'o', 'DisplayName','Homog C_{1122}')
if exist('Interpolation','var'); fplot(Interpolation.fun(1,1,2,2), [0 1], 'DisplayName','Fit/Interp'); end
fplot(C1122_simp, [0 1], 'DisplayName','SIMP')
xlabel('\phi'); ylabel('C_{1122}'); legend


nexttile; hold on; grid on
plot(volFrac, squeeze(Chomog(1,2,1,2,:)), 'o', 'DisplayName','Homog  C_{1212}')
if exist('Interpolation','var'); fplot(Interpolation.fun(1,2,1,2), [0 1], 'DisplayName','Fit/Interp'); end
fplot(C1212_simp, [0 1], 'DisplayName','SIMP')
xlabel('\phi'); ylabel('C_{1212}'); legend







% load('HomogenizationResults4.mat');
% simp = @(x) x.^3;
% 
% tiledlayout(1,3)
% nexttile
% hold on
% plot(volFrac,squeeze(Chomog(1,1,1,1,:)),'LineStyle','none','Marker','o')
% fplot(Interpolation.fun(1,1,1,1),[0 1])
% fplot(simp,[0 1])
% 
% nexttile
% hold on
% plot(volFrac,squeeze(Chomog(1,1,2,2,:)),'LineStyle','none','Marker','o')
% fplot(Interpolation.fun(1,1,2,2),[0 1])
% fplot(simp,[0 1])
% 
% nexttile
% hold on
% plot(volFrac,squeeze(Chomog(1,2,1,2,:)),'LineStyle','none','Marker','o')
% fplot(Interpolation.fun(1,2,1,2),[0 1])
% fplot(simp,[0 1])