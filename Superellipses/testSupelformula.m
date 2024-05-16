% Test polar coordinates
theta = 0:0.01:2*pi;
a = 0.8;
b = 0.6;
n = 3;
r_parameter = 0.5:0.1:1;

% Polar coordinates
% r = r_parameter*(abs(cos(theta)/a).^n+abs(sin(theta)/b).^n).^(-1/n);
% figure
% polarplot(theta,r)

figure
for i = 1:length(r_parameter)
    % Cartesian coordinates
    rc = r_parameter(i)*(abs(cos(theta)/a).^n+abs(sin(theta)/b).^n).^(-1/n);
    x = rc.*cos(theta);
    y = rc.*sin(theta);
    
    plot(x,y,'k','LineWidth', 2);
    grid on;
    axis equal;
    hold on
end
xlim([-1 1])
ylim([-1 1])
ta = annotation('textarrow', [0.515 0.75], [0.515 0.75]);
ta.Color = [0 0.5 0.5];
ta.String = 'r';  
ta.FontSize = 15;