load('HomogenizationResultsReinforcedHexagon.mat');   

E0   = 1.0;
nu   = 0.30;
p    = 3;
Emin = 1e-6;

Efun = @(phi) Emin + (E0 - Emin).*phi.^p;

% ---- PLANE STRESS ----
C1111_simp = @(phi) Efun(phi)./(1 - nu^2);
C2222_simp = @(phi) Efun(phi)./(1 - nu^2);
C1122_simp = @(phi) Efun(phi).*nu./(1 - nu^2);
tiledlayout(1,3)

% C1111
nexttile; hold on; grid on
plot(volFrac, squeeze(Chomog(1,1,1,1,:)), 'o', 'DisplayName','Homog C_{1111}')
if exist('Interpolation','var')
    fplot(Interpolation.fun(1,1,1,1), [0 1], 'DisplayName','Fit/Interp');
end
fplot(C1111_simp, [0 1], 'DisplayName','SIMP');
xlabel('\phi'); ylabel('C_{1111}'); legend

% C2222
nexttile; hold on; grid on
plot(volFrac, squeeze(Chomog(2,2,2,2,:)), 'o', 'DisplayName','Homog C_{2222}')
if exist('Interpolation','var')
    fplot(Interpolation.fun(2,2,2,2), [0 1], 'DisplayName','Fit/Interp');
end
fplot(C2222_simp, [0 1], 'DisplayName','SIMP');
xlabel('\phi'); ylabel('C_{2222}'); legend

% C1122
nexttile; hold on; grid on
plot(volFrac, squeeze(Chomog(1,1,2,2,:)), 'o', 'DisplayName','Homog C_{1122}')
if exist('Interpolation','var')
    fplot(Interpolation.fun(1,1,2,2), [0 1], 'DisplayName','Fit/Interp');
end
fplot(C1122_simp, [0 1], 'DisplayName','SIMP');
xlabel('\phi'); ylabel('C_{1122}'); legend





%%%%SIMP%%%%%

% 
% % --- SIMP parameters ---
% E0   = 1.0;
% nu   = 0.30;
% p    = 3;
% Emin = 1e-6;
% 
% % SIMP Young's modulus
% Efun = @(phi) Emin + (E0 - Emin).*phi.^p;
% 
% % ---- Plane stress stiffnesses ----
% C1111_simp = @(phi) Efun(phi)./(1 - nu^2);
% C2222_simp = @(phi) Efun(phi)./(1 - nu^2);
% C1122_simp = @(phi) Efun(phi).*nu./(1 - nu^2);
% 
% % Interval for plotting
% phi_range = [0 1];
% 
% % --------- PLOTS ---------
% tiledlayout(1,3)
% 
% % C1111
% nexttile; hold on; grid on
% fplot(C1111_simp, phi_range, 'LineWidth', 2)
% xlabel('\phi'); ylabel('C_{1111}')
% title('SIMP  C_{1111}')
% 
% % C2222
% nexttile; hold on; grid on
% fplot(C2222_simp, phi_range, 'LineWidth', 2)
% xlabel('\phi'); ylabel('C_{2222}')
% title('SIMP  C_{2222}')
% 
% % C1122
% nexttile; hold on; grid on
% fplot(C1122_simp, phi_range, 'LineWidth', 2)
% xlabel('\phi'); ylabel('C_{1122}')
% title('SIMP  C_{1122}')
