% Test polar coordinates
theta = 0:0.01:2*pi;
a = 1;
b = 1;
m = 5;
n1 = 2;
n2 = 6;
n3 = 6;

% Polar coordinates
r = (abs(cos(m*theta/4)/a).^n2+abs(sin(m*theta/4)/b).^n3).^(-1/n1);
figure
polarplot(theta,r)

% Cartesian coordinates
rc = (abs(cos(m*theta/4)/a).^n2+abs(sin(m*theta/4)/b).^n3).^(-1/n1);
x = rc.*cos(theta);
y = rc.*sin(theta);

figure
plot(x,y);
grid on;