function LevelSetVidergauz

E1 = 1e-8;
E2 = 1;
nu1 = 1/3;
nu2 = 1/3;

k = @(E,nu) E/(2*(1-nu));
mu = @(E,nu) E/(2*(1+nu));


k1 = k(E1,nu1);
k2 = k(E2,nu2);
mu1 = mu(E1,nu1);
mu2 = mu(E2,nu2);


[x,y] = createGrid();

theta = 0.2;
phi = 80*pi/180;

txi(1) = cos(phi);
txi(2) = sin(phi);
txi(3) = 0;

cx = 0.5;
cy = 0.5;


b = computeB(txi,mu2);
c = computeC(theta,txi,k1,k2,mu2);
q = b/c;

q2 = ((mu2 + k1)/(k2 - k1) + theta)*(txi(1) - txi(2))/(txi(1)+txi(2));

ax = computeAx(theta,q,cx,cy);
ay = computeAy(theta,q,cx,cy);

s.ax = ax;
s.ay = ay;
s.cx = cx;
s.cy = cy;

p = VigdergauzParametersComputer(s);

xp(:,1) = (ellipj(x(:,1)/p.rx*p.FxMax,p.mx));
yp(:,1) = (ellipj(y(:,1)/p.ry*p.FyMax,p.my));
levelset(:,1) = (1-xp(:,1).^2).*(1-yp(:,1).^2) - p.M;

hold on
isPositive = levelset>0;
validx = abs(x) <= p.rx;
validy = abs(y) <= p.ry;
ind = isPositive & validx & validy;
plot(x(~ind,:),y(~ind,:),'+')
axis([-cx, cx, -cy, cy])
end

function [x,y] = createGrid()
n = 1000;
x(:,1) = linspace(-0.5,0.5,n);
y(:,1) = linspace(-0.5,0.5,n);

x = repmat(x,1,500);
y = repmat(y',500,1);

x = x(:);
y = y(:);
end

function B = computeB(txi,mu2)
B = mu2*(txi(2) - txi(1) + 2*1i*txi(3));
end

function C = computeC(theta,txi,k1,k2,mu2)
C = (mu2*(k1 - k2)*(txi(1) + txi(2)))/(mu2+(1-theta)*k1+theta*k2);
end

function ax = computeAx(theta,q,cx,cy)
ax = (1 + theta + q)/4;
end

function ay = computeAy(theta,q,cx,cy)
ay = (1 + theta - q)/4;
end
