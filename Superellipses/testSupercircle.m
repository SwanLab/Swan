%% Program to obtain different type of supercircles

% Generate data points
x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);

a = 0.4;
b = 0.2;
p = 2;



    lhs = getSupercircleData(a,b,p,X,Y);
    contour(X, Y, lhs, [1 1], 'LineWidth', 2, 'LineColor', 'k');
    xlabel('x');
    ylabel('y');
    grid on;
    axis equal;
    title(['a = ', num2str(a),'   b = ', num2str(b),'   p = ', num2str(p)],'FontSize', 15)
    xlim([0 1])
    ylim([0 1])
