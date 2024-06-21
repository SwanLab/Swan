clc
clear
close all

%% Script to visualize the advection velocity
% Definition of the range of points
x_val = linspace(-0.5, 0.5, 20); 
y_val = linspace(-0.5, 0.5, 20); 

% Creation of the mesh (X, Y)
[X, Y] = meshgrid(x_val, y_val);

% Definition of advection velocity
% a = @(x, y) deal(-y ./ (x.^2 + y.^2), x ./ (x.^2 + y.^2));
% a = @(x, y) deal(-y, x);
% a = @(x,y) deal(-y+alpha.*x, x+alpha.*y);
% a = @(x,y) deal(alpha.*x, -alpha.*y);
% a = @(x,y) deal(-sin(x).*sin(y), cos(x).*cos(y));
% a = @(x,y) deal(2.*y.*exp(-(x.^2+y.^2)), -2.*x.*exp(-(x.^2+y.^2)));
% a = @(x,y) deal(sin(x).*cos(y), -cos(x).*sin(y));
% a = @(x,y) deal(x.^2-y.^2, -2.*x.*y);
a = @(x,y) deal(2*y.*cos(x.^2+y.^2), -2*x.*cos(x.^2+y.^2));

% Initialize u and v with the same size as X and Y
u = zeros(size(X));
v = zeros(size(Y));

% Calculate the velocity components at each point
for i = 1:numel(X)
    [u(i), v(i)] = a(X(i), Y(i));
end

% Plot the velocity field using quiver
figure;
q = quiver(X, Y, u, v);
q.Color = 'black';
title('Advection field');
xlabel('X');
ylabel('Y');
axis equal; % To ensure the scales of the axes are equal
grid on; % Activate the grid
xlim([-0.5 0.5])
ylim([-0.5 0.5])