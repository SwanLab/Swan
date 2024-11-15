function LaunchKinematics




g = 9.81;
y0 = 0;
dy0 = 100;

dx0 = 1;

Lf = 100;
error = 1;

nmax = 1e6;
while (itcount <= nmax && error >= 1e-6) 
itcount = itcount + 1; 
% Generate and save iteratres 
x = a + (b-a)/2; 
z(itcount) = x; 
fa = computeError(a,Lf);  
fb = computeDistance(b,Lf); 
fx = computeDistance(x,l);  
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




end
end

function error = computeError(dy0,Lf)
L = computeDistance(dy0);
error = L - Lf;
end

function [L] = computeDistance(dy0)
fun = @(t,y)  LaunchFunction(t,y);
tspan = [0 40];
yT0 = [y0,dy0];
[t,y] = ode45(fun, tspan, yT0);

posY = y(:,1)>0;

dx0 = 1;
L = dx0*t(posY);
%plot(t(posY),y(posY,1))


end

function LaunchFunction(t,y)
f = [y(2);-g];
end