%% Program to obtain different type of supercircles

% Generate data points
x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);

a = 1;
b = 0.5;
p = [0.5, 2, 5];
label = {'Hipoellipse','Ellipse','Hiperellipse'};
tiledlayout(1,3)

for i = 1:length(p)
    nexttile;
    lhs = getSupercircleData(a,b,p(i),X,Y);
    contour(X, Y, lhs, [1 1], 'LineWidth', 1.8, 'LineColor', 'k');
    xlabel('x');
ylabel('y');
grid on;
axis equal;
title([label{i}],'FontSize', 15)
xlim([-2 2])
ylim([-2 2])
end
