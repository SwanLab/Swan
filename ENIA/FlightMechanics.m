function FlightMechanics()


x0 = 0;
h0 = 10;
v0 = 2;
gamma0 = pi/8;
w0 = 1;

g = 9.81;
cP = 0.01;
Cla = 2*pi;
alpha = 0.1;

y0 = [x0;h0;v0;gamma0;w0];
tspan = [0:0.001:1];
fun = @(t,y) FlightMechanicsFunction(t,y,g,cP,Cla,alpha)

[t,y] = ode45(fun, tspan, y0);

x     = y(:,1);
h     = y(:,2);
v     = y(:,3);
gamma = y(:,4);
w     = y(:,5);

end


function f = FlightMechanicsFunction(t,y,g,cP,Cla,alpha)


h     = y(2);
v     = y(3);
gamma = y(4);
w     = y(5);
T = Thrust(h);
D = Drag(h,v,Cla,alpha);
L = Lift(h,v,Cla,alpha);

dx = v*cos(gamma);
dy = v*sin(gamma);
dv = g*((T-D)/w - sin(gamma));
dg = g/v*(L/w - cos(gamma));
dw = -cP*T;
f = [dx;dy;dv;dg;dw];
end

function T  = Thrust(h)
T = 0.7;
end

function L = Lift(h,v,Cla,alpha)
cL = Cla*alpha;
rho = Density(h);
S = 1;
L = 0.5*rho.*S*v.^2*cL;
end

function D = Drag(h,v,Cla,alpha)
cD = computeCd(Cla,alpha);
rho = Density(h);
S = 1;
D = 0.5*rho.*S*v.^2*cD;
end

function cD = computeCd(Cla,alpha)
Cd0 = 0.01;
eta = 0.95;
cD = Cd0+eta*Cla*alpha^2;
end

function rho = Density(h)
rho = 1;
end