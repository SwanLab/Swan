function R = test_b_value(obj, label)
%TEST_B_VALUE  Quantifies what the design variable b actually buys.
%
%   R = test_b_value(ans, 'a = 1')
%
%   Three questions, all answered from the database alone:
%
%   (1) GAIN        how much stiffness does b buy, per component, relative to
%                   the b = 0 cell? Measured against Hashin-Shtrikman, which is
%                   the attainable bound. Ratios to rho*C_solid (the straight
%                   line between void and solid) are NOT used: that line is not
%                   attainable and such ratios only restate the convexity of C.
%
%   (2) SATURATION  is the optimal b inside the admissible box, or is the box
%                   cutting off the optimum? If b_opt sits on the boundary for
%                   most of the rho range, widening the box still pays.
%
%   (3) CONFLICT    bulk and shear generally want opposite values of b. The
%                   angle between the two optima is the reason b has to be a
%                   spatial field rather than one global number: each point
%                   picks according to its local stress state.

if nargin < 2, label = ''; end
if nargin < 1 || isempty(obj), obj = evalin('base','ans'); end

Chomog   = obj.Chomog;
paramB   = obj.paramB(:).';
paramRho = obj.paramRho(:).';
nB = numel(paramB);  nRho = numel(paramRho);

E = 1; nu = 0.3;
k1  = E/(2*(1-nu));          % 2D bulk modulus of the solid
mu1 = E/(2*(1+nu));          % shear modulus of the solid

% Hashin-Shtrikman upper bounds in 2D (solid + void)
kHS  = @(r) r.*k1.*mu1 ./ ((1-r).*k1 + mu1);
muHS = @(r) r.*mu1.*(k1+2*mu1) ./ ((k1+2*mu1) + 2*(1-r).*(k1+mu1));

[~, ib0] = min(abs(paramB));      % the b closest to zero

fprintf('\n============================================================\n');
if ~isempty(label), fprintf('  %s\n', label); end
fprintf('  b in [%.2f, %.2f],  reference cell at b = %.4f\n', ...
        paramB(1), paramB(end), paramB(ib0));
fprintf('============================================================\n');

%% ------------------------------------- gather K and mu over the whole grid
K  = zeros(nRho, nB);
MU = zeros(nRho, nB);
for ir = 1:nRho
    for ib = 1:nB
        C = squeeze(Chomog(:,:,:,:,ir,ib));
        K(ir,ib)  = (C(1,1,1,1) + C(2,2,2,2) + 2*C(1,1,2,2))/4;
        MU(ir,ib) = C(1,2,1,2);
    end
end

%% ----------------------------------------------------------- (1) GAIN
fprintf('\n(1) what b buys, and how close it gets to the attainable bound\n');
fprintf('    %6s | %8s %8s %7s | %8s %8s %7s\n', ...
        'rho','mu(b=0)','mu(best)','gain','effHS b0','effHS bst','gain');
sel = round(linspace(max(2,round(0.1*nRho)), round(0.9*nRho), 8));
R.rho = paramRho(sel);
R.gainMu = zeros(size(sel));  R.gainK = zeros(size(sel));
for k = 1:numel(sel)
    ir = sel(k);  r = paramRho(ir);
    [muBest, ibM] = max(MU(ir,:));
    mu0 = MU(ir,ib0);
    R.gainMu(k) = muBest/max(mu0,eps);
    fprintf('    %6.3f | %8.4f %8.4f %7.2f | %8.2f %8.2f %7.2f\n', ...
            r, mu0, muBest, muBest/max(mu0,eps), ...
            mu0/muHS(r), muBest/muHS(r), (muBest-mu0)/muHS(r));
    R.bOptMu(k) = paramB(ibM);
end

fprintf('\n    same for the bulk modulus\n');
fprintf('    %6s | %8s %8s %7s | %8s %8s\n', ...
        'rho','K(b=0)','K(best)','gain','effHS b0','effHS bst');
for k = 1:numel(sel)
    ir = sel(k);  r = paramRho(ir);
    [kBest, ibK] = max(K(ir,:));
    k0 = K(ir,ib0);
    R.gainK(k) = kBest/max(k0,eps);
    R.bOptK(k) = paramB(ibK);
    fprintf('    %6.3f | %8.4f %8.4f %7.2f | %8.2f %8.2f\n', ...
            r, k0, kBest, kBest/max(k0,eps), k0/kHS(r), kBest/kHS(r));
end

%% ------------------------------------------------------ (2) SATURATION
fprintf('\n(2) is the box cutting off the optimum?\n');
onEdge = 0;  nUse = 0;
for ir = 2:nRho-1
    if paramRho(ir) < 0.1 || paramRho(ir) > 0.9, continue; end
    nUse = nUse + 1;
    [~, ibM] = max(MU(ir,:));
    if ibM <= 2 || ibM >= nB-1, onEdge = onEdge + 1; end
end
fprintf('    optimal b for shear on the box boundary: %.0f%% of the rho range\n', ...
        100*onEdge/max(nUse,1));
R.edgeFrac = onEdge/max(nUse,1);
if onEdge/max(nUse,1) > 0.5
    fprintf('    -> the box is the binding constraint, not the physics.\n');
    fprintf('       Widen b and rerun: the shear stiffness has not peaked yet.\n');
else
    fprintf('    -> the optimum is interior: the admissible range is wide enough.\n');
end

%% -------------------------------------------------------- (3) CONFLICT
fprintf('\n(3) bulk and shear want different b\n');
fprintf('    %6s %10s %10s %10s\n','rho','b* shear','b* bulk','separation');
for k = 1:numel(sel)
    fprintf('    %6.3f %10.3f %10.3f %10.3f\n', ...
            R.rho(k), R.bOptMu(k), R.bOptK(k), abs(R.bOptMu(k)-R.bOptK(k)));
end
fprintf('    -> a nonzero separation is the argument for b as a SPATIAL FIELD:\n');
fprintf('       shear-dominated and compression-dominated regions of the same\n');
fprintf('       structure ask for different microstructures.\n');

%% -------------------------------------------------------------- figures
figure('Name',['value of b   ' label],'Position',[60 60 1200 420]);

subplot(1,3,1)
plot(paramRho, MU(:,ib0), 'k-', 'LineWidth',1.6); hold on
plot(paramRho, max(MU,[],2), 'r-', 'LineWidth',1.6);
plot(paramRho, muHS(paramRho), 'b--', 'LineWidth',1.2);
grid on; xlabel('\rho'); ylabel('shear modulus')
legend('b = 0','best b','HS bound','Location','northwest')
title('what b buys in shear')

subplot(1,3,2)
plot(paramRho, K(:,ib0), 'k-', 'LineWidth',1.6); hold on
plot(paramRho, max(K,[],2), 'r-', 'LineWidth',1.6);
plot(paramRho, kHS(paramRho), 'b--', 'LineWidth',1.2);
grid on; xlabel('\rho'); ylabel('bulk modulus')
legend('b = 0','best b','HS bound','Location','northwest')
title('what b buys in bulk')

subplot(1,3,3)
bOptMu = zeros(nRho,1);  bOptK = zeros(nRho,1);
for ir = 1:nRho
    [~,i1] = max(MU(ir,:));  bOptMu(ir) = paramB(i1);
    [~,i2] = max(K(ir,:));   bOptK(ir)  = paramB(i2);
end
plot(paramRho, bOptMu, 'r-', paramRho, bOptK, 'b-', 'LineWidth',1.6); grid on
yline(paramB(1),'k:'); yline(paramB(end),'k:');
xlabel('\rho'); ylabel('optimal b'); ylim([paramB(1) paramB(end)])
legend('shear','bulk','Location','best')
title('the two optima disagree')

R.K = K;  R.MU = MU;  R.paramB = paramB;  R.paramRho = paramRho;
R.label = label;
end