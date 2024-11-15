function Pendulum()


g = 9.81;
l = 3;
w = sqrt(g/l);




tspan = [0 50];
y0 = [pi/8;0];

fun2Integrate = @(t,y) [y(2);-w^2*sin(y(1))];

[t1,y1] = integrateOde(fun2Integrate, tspan, y0);
plotSolution(t1,y1,g,l)


[t1,y1] = ode45(fun2Integrate, tspan, y0);
plotSolution(t1,y1,g,l)


fun2Integrate = @(t,y) [y(2);-w^2*(y(1))];
[t2,y2] = ode45(fun2Integrate, tspan, y0);
plotSolution(t2,y2,g,l)

end

function dydtLinearized(obj,)

end

function dydtLinearized()

end



function [tv,y] = integrateOde(fun,tspan,y0)
h = 0.001;
tv = linspace(tspan(1),tspan(2),10000);
y(1,:) = y0;
for it = 1:length(tv)-1
    fv = fun(tv(it),y(it,:));
    y(it+1,:) = y(it,:) + fv';
end

end

function plotSolution(t,y,g,l)
figure
plot(t,y(:,1))

th = y(:,1);
dth= y(:,2);


Ep = (1 - cos(th))*g*l;
Ec = 0.5*dth.^2*l^2;
figure
plot(t,[Ep,Ec])

figure
plot(t,[Ep+Ec])

end
