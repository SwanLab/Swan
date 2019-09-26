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

mx = computeEllipticParameter(ax,ay,cy);
my = computeEllipticParameter(ay,ax,cx);

Kx = completeEllipticFunction(mx);
Ky = completeEllipticFunction(my);
M = computeM(mx,my);

FxMax = incompleteEllipticFunction(sqrt(1-M),mx);
FyMax = incompleteEllipticFunction(sqrt(1-M),my);

rx = ax/Kx*FxMax;
ry = ay/Ky*FyMax;

xp(:,1) = (ellipj(x(:,1)/rx*FxMax,mx));
yp(:,1) = (ellipj(y(:,1)/ry*FyMax,my));
levelset(:,1) = (1-xp(:,1).^2).*(1-yp(:,1).^2) - M;

hold on
isPositive = levelset>0;
validx = abs(x) <= rx;
validy = abs(y) <= ry;
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

function Tx = computeT(a1,a2,c)
%Tx = (1-theta + q)/(1+theta + q);
Tx = (c - a2)/a1;
end
% 
% function Ty = computeT(a1,a2,c)
% %Ty = (1-theta - q)/(1+theta - q);
% Ty = (c - a2)/a1;
% end

function K = completeEllipticFunction(k)
K = incompleteEllipticFunction(1,k);
end

function F = incompleteEllipticFunction(x,k)
F = ellipticF(asin(x),k);
x2 = (ellipj(F,k));
%norm(x -x2)
end

function x = computeEllipticParameter(a1,a2,c)
T = computeT(a1,a2,c);
F = @(x) implicitMequation(x,T);
eps = 1e-15;
x = fzero(F,[eps,1-eps]);
end

function f = implicitMequation(x,T)
f = completeEllipticFunction(x)*T - completeEllipticFunction(1-x);
end

function M = computeM(mx,my)
f = @(x) (1 - x)/x;
M = f(mx)*f(my);
end

