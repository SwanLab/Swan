function LaunchKinematics






dy0A = 100;
dy0B = 1;


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

dy0 = a;
fun = @(t,y)  LaunchFunction(t,y);
tspan = [0 40];
y0 = 0;
yT0 = [y0,dy0];
[t,y] = ode45(fun, tspan, yT0);
posY = y(:,1)>0;
plot(t(posY),y(posY,1))

end
end

function error = computeError(dy0,Lf)
L = computeDistance(dy0);
error = L - Lf;
end

function [L] = computeDistance(dy0)
y0 = 0;
fun = @(t,y)  LaunchFunction(t,y);
tspan = [0 40];
yT0 = [y0,dy0];
[t,y] = ode45(fun, tspan, yT0);

posY = y(:,1)>0;

dx0 = 1;
tpos = t(posY);
L = dx0*tpos(end);
%plot(t(posY),y(posY,1))


end

function f = LaunchFunction(t,y)
g = 9.81;
f = [y(2);-g];
end