function Lecture2Exercice2
t0   = 0;
tmax = 1;
nT = 100;
h = (tmax - t0)/nT;
t = linspace(0,1,nT);

g = 9.81;
L = 1;

theta0  = pi/20;
dtheta0 = 0;

y0 = [theta0; dtheta0];
y = zeros(2,nT);
y(:,1) = y0;
dy = @(t,y) computeDerivative(t,y,g,L);

for it = 1:nT-1
    yt1 = forwardEuler(t(it),y(:,it),h,dy);   
    y(:,it+1) = yt1;
end
plot(t,y)


end

function dy = computeDerivative(t,y,g,L)
dy(1,1) = y(2);
dy(2,1) = -g/L*sin(y(1));
end

function dy = computeLinearizedDerivative(t,y,g,L)
dy(1,1) = y(2);
dy(2,1) = -g/L*(y(1));
end

function yt1 = forwardEuler(t,yt,h,dy)
f = dy(t,yt);
yt1 = yt + h*f;
end