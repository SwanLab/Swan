function phi = LevelSetLipung(x,y,theta,phi)

x = (1 - (-1))*(x-0.5);
y = (1 - (-1))*(y-0.5);

cx = 1; cy = 1;



r = computeOptimalR(theta,phi,cx,cy);

h = computeH(theta,r);

Tx = computeTx(r,h,cx);
Ty = computeTy(r,h,cx);

ax = computeAx(Tx,Ty,cx);
ay = computeAy(Tx,Ty,cx);

mx = computeEllipticParameter(ax,ay,cy);
my = computeEllipticParameter(ay,ax,cx);

M = computeM(mx,my);

FxMax = incompleteEllipticFunction(sqrt(1-M),mx);
FyMax = incompleteEllipticFunction(sqrt(1-M),my);

rx = computeR(ax,mx,M);
ry = computeR(ay,my,M);

%[x,y] = createGrid(cx,cy);

xp(:,1) = (ellipj(x(:,1)/rx*FxMax,mx));
yp(:,1) = (ellipj(y(:,1)/ry*FyMax,my));
levelset(:,1) = (1-xp(:,1).^2).*(1-yp(:,1).^2) - M;

hold on
isPositive = levelset>0;
validx = abs(x) <= rx;
validy = abs(y) <= ry;
valid =  validx & validy;
ind = isPositive & valid;
%ind = isPositive;
%plot(x(~ind,:),y(~ind,:),'+')
%axis([-cx, cx, -cy, cy])

phi = levelset;

%phi(ind) = -phi(ind);
phi(~(validx & validy)) = -abs(phi(~(validx & validy)));

% tri = delaunay(x,y);
% h = trisurf(tri, x, y, phi);
% axis off
% shading interp

end

function Tx = computeTx(r,h,cx)
n = r*h*(1+h);
d = (1+r*h)*cx;
Tx = n/d;
end

function Ty = computeTy(r,h,cx)
n = cx*h*(1+r*h);
d = (1+h);
Ty = n/d;
end

function ax = computeAx(Tx,Ty,cx)
ax = (cx - Ty)/(1-Tx*Ty);
end

function ay = computeAy(Tx,Ty,cx)
ay = (1-Tx*cx)/(1-Tx*Ty);
end

function r = computeOptimalR(theta,phi,cx,cy)
F = @(r) equationForR(r,theta,phi,cx,cy);
eps = 0.5;
c = 1.1;
% ru = c;
% rl = 1/c;
[rub,rlb] = findRbounds2(F);
options = optimset('Display','iter');
%r = fzero(F,[rlb,rub],options);
r = fzero(F,[rlb,rub]);
%options = optimset('Display','iter');
%r = fzero(F,1,options);
end

function [rub,rlb] = findRbounds(F)
F1 = F(1);
if F1 >= 0   
    Fr = F1;
    i = -6;
    while Fr >= 0
       r = 1 + 10^(i);
       Fr = F(r);      
       i = i + 1;
    end
    rub = r;
    rlb = 1 + 10^(i-2);
else
    Fr = F1;
    i = -10;
    while Fr <= 0
       r = 1 - 10^(i);
       Fr = F(r);      
       i = i + 1;
    end
    rub = r;
    rlb = 1 - 10^(i-2);    
end


end

function [rub,rlb] = findRbounds2(F)
r0 = 1;
F0 = F(1);
n = 6;
if F0 >= 0   
    r1 = 1 + 10^(-n);
    F1 = F(r1);   
    while F1 >= 0
       rnew = newPointBySecant(r0,r1,F0,F1);
       r0 = r1;
       F0 = F1;
       r1 = rnew;
       F1 = F(r1);      
    end
    rub = r1;
    rlb = r0;    
else
    r0 = 1;
    r1 = 1 - 10^(-n);
    F1 = F(r1);   
    while F1 <= 0
       rnew = newPointBySecant(r0,r1,F0,F1);
       r0 = r1;
       F0 = F1;
       r1 = max(1e-12,rnew);
       %r1 = rnew;
       F1 = F(r1);      
    end
    rub = r0;
    rlb = r1;        
end

end

function x2 = newPointBySecant(x0,x1,f0,f1)
x2 = x1 - (x1-x0)/(f1 - f0)*f1;
end


function f = equationForR(r,theta,phi,cx,cy)

h = computeH(theta,r);

Tx = computeTx(r,h,cx);
Ty = computeTy(r,h,cx);

ax = computeAx(Tx,Ty,cx);
ay = computeAy(Tx,Ty,cx);

mx = computeEllipticParameter(ax,ay,cy);
my = computeEllipticParameter(ay,ax,cx);

M = computeM(mx,my);
rx = computeR(ax,mx,M);
ry = computeR(ay,my,M);

f = tan(phi) - rx/ry;
end

function h = computeH(theta,r)
n = -(1+r)*theta + sqrt(theta^2*(r-1)^2 + 4*r);
d = 2*r*(1+theta);
h = n/d;
end

function r = computeR(a,m,M)
K = completeEllipticFunction(m);
Fmax = incompleteEllipticFunction(sqrt(1-M),m);
r = a/K*Fmax;
end

function x = computeEllipticParameter(a1,a2,c)
T = computeT(a1,a2,c);
eps = 2.2204*1e-14;
if T <= 0.15
    x = 1 - eps;
else
    F = @(x) implicitMequation(x,T);
    options = optimset('Display','iter');
    %x = fzero(F,[0+eps,1-eps],options);
    x = fzero(F,[0+eps,1-eps]);
end
end

function f = implicitMequation(x,T)
f = completeEllipticFunction(x)*T - completeEllipticFunction(1-x);
end

function K = completeEllipticFunction(k)
K = incompleteEllipticFunction(1,k);
end

function F = incompleteEllipticFunction(x,k)
F = ellipticF(asin(x),k);
end

function Tx = computeT(a1,a2,c)
%Tx = (1-theta + q)/(1+theta + q);
Tx = (c - a2)/a1;
end

function M = computeM(mx,my)
f = @(x) (1 - x)/x;
M = f(mx)*f(my);
end


function [x,y] = createGrid(cx,cy)
n = 1000;
x(:,1) = linspace(-cx,cx,n);
y(:,1) = linspace(-cy,cy,n);

x = repmat(x,1,500);
y = repmat(y',500,1);

x = x(:);
y = y(:);
end