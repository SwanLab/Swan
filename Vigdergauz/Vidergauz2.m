function Vidergauz2

myCase = 'B';

switch myCase
    case 'A'
        nTheta = 10;
        ntheta0 = -3;
        nthetaF = -0.1;
        %theta1V = linspace(theta0,thetaF,nTheta);
        thetaV = 10.^linspace(ntheta0,nthetaF,nTheta);
        
        rV = 1;
        cxV = 1;
        cy = 1;
        ncx = length(cxV);
    case 'B'
        nTheta = 1;
        thetaV = 0.5;
        
        r0 = linspace(80,2,5);
        rf = 1./r0;
        rV = sort([r0 1 rf]);
       % rV = [1/5 1/4 1/3 1/2 1 2 3 4 5];
        
        cxV = 1;
        cy = 1;
        ncx = length(cxV);
        
    case 'C'
        thetaV = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
        %thetaV = 0.5/4;
        nTheta = length(thetaV);
        rV = 1;
        cxV = 1;
        cy = 1;
        ncx = length(cxV);
        
    case 'D'
        nTheta = 1;
        thetaV = 0.5;
        
        rV = 1;
        cxV = 1;
        cy = 1;
        ncx = length(cxV);        
        
        
end



figure(1)
hold on
for itheta = 1:nTheta
    theta = thetaV(itheta);
    for icx = 1:ncx
        cx = cxV(icx);
        for ir = 1:length(rV)
         r = rV(ir);
        
        h = computeH(theta,r);
        
        ax = computeAx(h,cx)
        ay = computeAy(h,r)
        
        mx = computeSmallM(ax,ay,cy);
        my = computeSmallM(ay,ax,cx);
        
        ax2 = (2*r*(theta + 1))/(2*r - theta + (r^2*theta^2 - 2*r*theta^2 + 4*r + theta^2)^(1/2) + r*theta);
        ay2 = (2*theta + 2)/(theta + (r^2*theta^2 - 2*r*theta^2 + 4*r + theta^2)^(1/2) - r*theta + 2);

        ax3 = (cx*cy - cx + cx*r + cx*cy*theta)/(cy + r);
        ay3 = (cy + cy*r*theta)/(cy + r);
        
        Tx = computeTx(r,h,cx);
        Ty = computeTy(r,h,cx);
        
        ax4 = (cx - Ty)/(1-Tx*Ty);
        ay4 = (1 - Tx*cx)/(1-Tx*Ty);
        
        mx2 = computeSmallMx(ax,ay,cy);
        my2 = computeSmallMy(ax,ay,cx);
        
        h2 = (1-ax/cx)/ay;
        r2 = (1-ay)/(1-ax/cx);
        
        r3 = (ax - ax*ay)/(ay - ax*ay);
        
        M = computeM(mx,my);
        
        Kx = ellipticF(asin(1),(mx));
        Ky = ellipticF(asin(1),(my));
        
        FxMax = incompleteEllipticFunction(sqrt(1-M),mx);
        FyMax = incompleteEllipticFunction(sqrt(1-M),my);
        
        rx = ax/Kx*FxMax;
        ry = ay/Ky*FyMax;
        
        [xq,yq] = computeShapePoints(ax,ay,mx,my,Kx,Ky,M);
        
        plotShape(xq,yq);
        end
    end
end

end

function Tx = computeT(a1,a2,c)
Tx = (c - a2)/a1;
end

function my = computeSmallMy(ax,ay,cx)
f = @(x) cx - ax - ay*ellipticF(asin(1),(1-x))/ellipticF(asin(1),x);
eps = 1e-14;
my = fzero(f,[eps,1-eps]);
end

function mx = computeSmallMx(ax,ay,cy)
f = @(x) cy - ax*ellipticF(asin(1),(1-x))/ellipticF(asin(1),x) - ay;
eps = 1e-14;
mx = fzero(f,[eps,1-eps]);
end

function h = computeH(theta,r)
n = -(1+r)*theta + sqrt(theta^2*(r-1)^2 + 4*r);
d = 2*r*(1+theta);
h = n/d;
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

function ax = computeAx(h,cx)
ax = cx/(1+h);
end

function ay = computeAy(h,r)
ay = 1/(1+r*h);
end

function m = computeSmallM(a1,a2,c)
T = computeT(a1,a2,c);
f = @(x) ellipticF(asin(1),x)*T - ellipticF(asin(1),(1-x));
eps = 1e-14;
m = fzero(f,[eps,1-eps]);
end

function M = computeM(mx,my)
n = (1-mx)*(1-my);
d = mx*my;
M = n/d;
end

function x = computeX(t,ax,mx,Kx)
xe = sqrt(1-t);
F = ellipticF(asin(xe),mx);
x = ax*F/Kx;
%x = ax*F/2;
end

function y = computeY(t,ay,my,Ky,M)
ye = sqrt(1-M/t);
F = ellipticF(asin(ye),my);
y = ay*F/Ky;
%y = ay*F/2;
end

function [xq,yq] = computeShapePoints(ax,ay,mx,my,Kmx,Kmy,M)
Nt = 100;
t = linspace(M,1,Nt);
for it = 1:Nt
    x(it,1) = computeX(t(it),ax,mx,Kmx);
    y(it,1) = computeY(t(it),ay,my,Kmy,M);
end
xq = [x;x;-x;-x];
yq = [-y;y;y;-y];
end

function plotShape(xq,yq)
plot(xq,yq,'+','MarkerSize',20)
h = 1;
axis([-h, h, -h, h])
drawnow
end

function F = incompleteEllipticFunction(x,k)
F = ellipticF(asin(x),k);
x2 = (ellipj(F,k));
%norm(x -x2)
end