% Test polar coordinates
theta = 0:0.01:2*pi;
a = 0.8;
b = 0.6;
p = 0.75;

% Polar coordinates
% r = r_parameter*(abs(cos(theta)/a).^n+abs(sin(theta)/b).^n).^(-1/n);
% figure
% polarplot(theta,r)

r = (abs(cos(theta)).^p+abs(sin(theta)).^p).^(-1/p);

x = r.*cos(theta);
y = r.*sin(theta);

t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(theta,r,'k','LineWidth', 2);
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlim([0 2*pi])
ylim([0 3])
xlabel('\theta')
ylabel('\rho')
grid on;
nexttile
polarplot(theta,r,'k','LineWidth', 2)
thetaticks([0 30 60 90 120 150 180 210 240 270 300 330])
thetaticklabels({'0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6','\pi','7\pi/6','4\pi/3','3\pi/2','5\pi/3','11\pi/6'})
title(t,['p = ', num2str(p), '    ', 'f_{min} = ', num2str(min(r)), '    ', 'f_{max} = ', num2str(max(r))])