%% Program to obtain different type of supercircles

% Generate data points
x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);

a = 1;
b = 0.5;
p = [0.5, 2, 10];
label = {'Hipoellipse', 'Ellipse', 'Superellipse'};

t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for i = 1:3
    nexttile
    lhs = getSupercircleData(a,b,p(i),X,Y);
    contour(X, Y, lhs, [1 1], 'LineWidth', 2, 'LineColor', 'k');
    xlabel('x');
    ylabel('y');
    title(label{i}, 'FontSize', 20);
    grid on;
    axis equal;
end



