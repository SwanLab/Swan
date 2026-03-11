%% TEST_RHOMBOID_ALPHA - usando TutorialHomogenization

clear; close all; clc;

c1 = 1.0;
Area = 1.0;
alpha_vals = 30:5:150;

nAlpha = length(alpha_vals);

C11_vals = zeros(nAlpha,1);
C22_vals = zeros(nAlpha,1);
C12_vals = zeros(nAlpha,1);
C33_vals = zeros(nAlpha,1);
vf_vals  = zeros(nAlpha,1);

for i = 1:nAlpha
    
    alpha = alpha_vals(i);
    c2 = Area / (c1*sind(alpha));
    
    fprintf('α = %4d° → c2 = %.4f ', alpha, c2);
    
    try
        TH = TutorialHomogenization(c1,c2,alpha,40);
        [C,vf] = TH.runAtVolume(0.5);
        
        C11_vals(i) = C(1,1,1,1);
        C22_vals(i) = C(2,2,2,2);
        C12_vals(i) = C(1,1,2,2);
        C33_vals(i) = C(1,2,1,2);
        vf_vals(i)  = vf;
        
        fprintf('→ vf=%.3f\n', vf);
        
    catch ME
        fprintf('→ ERROR: %s\n', ME.message);
        C11_vals(i) = NaN;
        C22_vals(i) = NaN;
        C12_vals(i) = NaN;
        C33_vals(i) = NaN;
        vf_vals(i)  = NaN;
    end
end

% 🔥 EXATAMENTE o mesmo plot que você já tinha
plotAlphaResults(alpha_vals, ...
                 C11_vals, C22_vals, ...
                 C12_vals, C33_vals, ...
                 vf_vals);
%% =====================================================
%% FUNÇÃO DE PLOT
%% =====================================================

function plotAlphaResults(alpha_vals, C11, C22, C12, C33, vf)

    colors = [
        0.2667, 0.4667, 0.6667;  % blue
        0.9333, 0.4000, 0.4667;  % red
        0.4000, 0.7333, 0.4000;  % green
        0.9333, 0.6667, 0.2000;  % orange
        0.6000, 0.4000, 0.6667;  % purple
    ];

    figure('Color','white','Position',[100 100 1400 1000])

    subplot(2,2,1)
    plot(alpha_vals,C11,'o-','LineWidth',1.5,'Color',colors(1,:))
    xlabel('\alpha (deg)')
    ylabel('C_{11}')
    grid on

    subplot(2,2,2)
    plot(alpha_vals,C22,'s-','LineWidth',1.5,'Color',colors(2,:))
    xlabel('\alpha (deg)')
    ylabel('C_{22}')
    grid on

    subplot(2,2,3)
    plot(alpha_vals,C33,'d-','LineWidth',1.5,'Color',colors(3,:))
    xlabel('\alpha (deg)')
    ylabel('C_{33}')
    grid on

    subplot(2,2,4)
    plot(alpha_vals,C12,'^-','LineWidth',1.5,'Color',colors(4,:))
    xlabel('\alpha (deg)')
    ylabel('C_{12}')
    grid on

    sgtitle('Homogenized Tensor vs Angle \alpha (vf ≈ 0.5)')

    % Anisotropy ratio
    figure('Color','white')
    plot(alpha_vals,C11./C22,'*-','LineWidth',2,'Color',colors(5,:))
    xlabel('\alpha (deg)')
    ylabel('C_{11}/C_{22}')
    grid on
    yline(1,'--','Isotropy')

    % Volume fraction check
    figure('Color','white')
    plot(alpha_vals,vf,'o-','LineWidth',1.5)
    xlabel('\alpha (deg)')
    ylabel('vf')
    grid on
    yline(0.5,'--','Target')
end