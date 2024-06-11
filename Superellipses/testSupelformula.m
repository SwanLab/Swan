% Test polar coordinates
n = 2;
theta = 0:0.01:n*2*pi;
a = 10;
b = 10;
p1 = 2;
p2 = 4.6;
p3 = 5.9;
m = 4;

% t = tiledlayout(1,5);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';

r = (abs(1/a*cos(m/4*theta)).^p2+abs(1/b*sin(m/4*theta)).^p3).^(-1/p1);

x = r.*cos(theta);
y = r.*sin(theta);

% nexttile
% plot(theta,r,'k','LineWidth', 2);
% xticks([0 pi/2 pi 3*pi/2 2*pi])
% xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
% xlim([0 2*pi])
% ylim([0 3])
% xlabel('\theta')
% ylabel('\rho')
% grid on;
nexttile
polarplot(theta,r,'k','LineWidth', 2)
thetaticks([0 30 60 90 120 150 180 210 240 270 300 330])
thetaticklabels({'0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6','\pi','7\pi/6','4\pi/3','3\pi/2','5\pi/3','11\pi/6'})
title({['p_1 = ', num2str(p1),'  p_2 = ', num2str(p2),'  p_3 = ', num2str(p3)], ['m = ', num2str(m), '  n = ', num2str(n), ' rev']}, 'FontSize', 15);


