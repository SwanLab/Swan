clear all
g = sym('g','positive');
l = sym('l','positive');

w = sqrt(g/l)

 A = sym('A');
 B = sym('B');
% 
 t = sym('t');
% 
 theta0 = sym('theta0','real');
 dtheta0 = sym('dtheta0','real');
% 
% 
% 
 theta = A*cos(w*t) + B*sin(w*t);
% 
dtheta = diff(theta,t)
% 
 eq(1) = theta0  == subs(theta,t,0)   
 eq(2) = dtheta0 == subs(dtheta,t,0)
% 
 solution = solve(eq,[A,B])
% 
 theta = subs(theta,[A,B],[solution.A,solution.B]);
 theta = simplify(theta);
% 

theta = subs(theta,g,9.81);
theta = subs(theta,l,1);
theta = subs(theta,theta0,pi/8);
theta = subs(theta,dtheta0,0)

theta = matlabFunction(theta);
tV = linspace(0,1,100);
theta(tV)
plot(tV,theta(tV),'+')

