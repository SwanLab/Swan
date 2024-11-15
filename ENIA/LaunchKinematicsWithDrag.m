function LaunchKinematicsWithDrag



dy0A = 1;
dy0B = 0.1;


[LA] = computeDistance(dy0A);
[LB] = computeDistance(dy0B);
Lf = 1;
error = 1;

a = dy0A;
b = dy0B;
nmax = 1e6;
itcount = 1;
tol =1e-6;
while (itcount <= nmax && error >= tol) 
itcount = itcount + 1; 
% Generate and save iteratres 
x = a + (b-a)/2; 
z(itcount) = x; 
fa = computeError(a,Lf);  
fb = computeError(b,Lf); 
fx = computeError(x,Lf);  
error = abs(fx); 
%  error = abs(x - xold); 
if (error < tol) 
x_final = x; 
else 
if (fa*fx < 0) 
% root is between a and x 
b = x; 
else 
% root is between x and b 
a = x; 
end 
end 

[Ls] = computeDistance(x);

v0 = a;
fun = @(t,y)  LaunchFunction(t,y);
tspan = [0:0.01:40];
x0 = 0;
y0 = 0;
theta0 = pi/8;
yT0 = [x0,y0,v0,theta0];
[t,y] = ode45(fun, tspan, yT0);
posY = y(:,2)>0;
plot(t(posY),y(posY,1))

end
end

function error = computeError(dy0,Lf)
L = computeDistance(dy0);
error = L - Lf;
end

function [L] = computeDistance(theta0)
x0 = 0;
y0 = 0;
v0 = 1;
fun = @(t,y)  LaunchFunction(t,y);
tspan = [0:0.01:40];
yT0 = [x0,y0,v0,theta0];
[t,y] = ode45(fun, tspan, yT0);

posY = y(:,2)>0;

dx0 = 1;
tpos = t(posY);
L = dx0*tpos(end);
%plot(t(posY),y(posY,1))


end

function f = LaunchFunction(t,y)
mu = 0;
xv = y(1);
yv = y(2);
v  = y(3);
theta = y(4);
g = 9.81;
f = [v*cos(theta);v*sin(theta);-g*sin(theta)-mu*v^2;-g*cos(theta)/v];
end