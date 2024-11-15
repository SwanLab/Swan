function Pendulum()
% clear all
% g = sym('g','positive');
% l = sym('l','positive');
% 
% w = sqrt(g/l)
% 
%  A = sym('A');
%  B = sym('B');
% % 
%  t = sym('t');
% % 
%  theta0 = sym('theta0','real');
%  dtheta0 = sym('dtheta0','real');
% % 
% % 
% % 
%  theta = A*cos(w*t) + B*sin(w*t);
% % 
% dtheta = diff(theta,t)
% % 
%  eq(1) = theta0  == subs(theta,t,0)   
%  eq(2) = dtheta0 == subs(dtheta,t,0)
% % 
%  solution = solve(eq,[A,B])
% % 
%  theta = subs(theta,[A,B],[solution.A,solution.B]);
%  theta = simplify(theta);
% % 
% 
% theta = subs(theta,g,9.81);
% theta = subs(theta,l,1);
% theta = subs(theta,theta0,pi/8);
% theta = subs(theta,dtheta0,0)
% 
% theta = matlabFunction(theta);
% tV = linspace(0,1,100);
% theta(tV)
% plot(tV,theta(tV),'+')

g = 9.81;
l = 3;
w = sqrt(g/l);




tspan = [0 50];
y0 = [pi/8;0];

odefun = @(t,y) [y(2);-w^2*sin(y(1))];

[t1,y1] = integrateOde(odefun, tspan, y0);
plotSolution(t1,y1,g,l)


[t1,y1] = ode45(odefun, tspan, y0);
plotSolution(t1,y1,g,l)


odefun = @(t,y) [y(2);-w^2*(y(1))];
[t2,y2] = ode45(odefun, tspan, y0);
plotSolution(t2,y2,g,l)

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
