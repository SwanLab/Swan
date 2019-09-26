function Vigdergauz
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


myCase = 'B';

switch myCase
    case 'A'
        nTheta = 5;
        ntheta0 = -3;
        nthetaF = -0.1;
        %theta1V = linspace(theta0,thetaF,nTheta);
        thetaV(:,1) = 10.^linspace(ntheta0,nthetaF,nTheta)
        
        nPhi = 1;
        phiT = 45;
    case 'B'
        nTheta = 1;
        thetaV = 0.1;
        
        nPhi = 5;
        phi0 = 0;
        phi1 = 90;
        phiT = linspace(phi0,phi1,nPhi);
        
    case 'C'
        thetaV = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
        nTheta = length(thetaV);
        
        nPhi = 1;
        phiT = 45;
end


figure(1)
hold on

cx = 1/2;
cy = 1/2; 

for itheta = 1:nTheta
    theta = thetaV(itheta);
    for iPhi = 1:nPhi
        
        phi = phiT(iPhi)*(pi/180);
        txi(1) = cos(phi);
        txi(2) = sin(phi);
        txi(3) = 0;
        
        b = computeB(txi,mu2);
        c = computeC(theta,txi,k1,k2,mu2);
        q = b/c;
        
        
        
        ax = computeAx(theta,q,cx,cy);
        ay = computeAy(theta,q,cx,cy);
        
        Tx = computeTx(theta,q);       
        Ty = computeTy(theta,q);
        
        Tx2 = (cy - ay)/ax;
        Ty2 = (cx - ax)/ay;
        
        
        isValid = Ty - 1/Tx < 0
        
        mx = computeEllipticParameter(Tx);
        my = computeEllipticParameter(Ty);
        
        M = computeM(mx,my);
        
        Kx = completeEllipticFunction(mx);
        Ky = completeEllipticFunction(my);
        
        FxMax = incompleteEllipticFunction(sqrt(1-M),mx);
        FyMax = incompleteEllipticFunction(sqrt(1-M),my);
        
        rx = ax/Kx*FxMax;
        ry = ay/Ky*FyMax;
        
        [ax rx; ay ry]
        
        [xq,yq] = computeVidgergauzCurve(mx,my,ax,ay);
        plotVigdergauz(xq,yq,cx,cy)
    end
end

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

function Tx = computeTx(theta,q)
Tx = (1-theta + q)/(1+theta + q);
end

function Ty = computeTy(theta,q)
Ty = (1-theta - q)/(1+theta - q);
end

function x = computeEllipticParameter(p)
F = @(x) implicitMequation(x,p);
eps = 1e-15;
x = fzero(F,[eps,1-eps]);
end

function [x,y] = computeVidgergauzCurve(mx,my,ax,ay)
Kx = completeEllipticFunction(mx);
Ky = completeEllipticFunction(my);
M = computeM(mx,my);
[xq,yq] = computeQuarterCurvePoints(ax,ay,mx,my,Kx,Ky,M);
[x,y] = computeAllCurvePoints(xq,yq);
end

function [x,y] = computeQuarterCurvePoints(ax,ay,mx,my,Kx,Ky,M)
Nt = 100;
t = linspace(M,1,Nt);
x = zeros(Nt,1);
y = zeros(Nt,1);
for it = 1:Nt
    x(it,1) = computeX(ax,Kx,mx,t(it));
    y(it,1) = computeY(ay,Ky,my,M,t(it));
end
end

function [x,y] = computeAllCurvePoints(xq,yq)
x = [xq;xq;-xq;-xq];
y = [-yq;yq;yq;-yq];
end

function M = computeM(mx,my)
f = @(x) (1 - x)/x;
M = f(mx)*f(my);
end

function plotVigdergauz(xq,yq,cx,cy)
plot(xq,yq,'+','MarkerSize',20)
axis([-cx, cx, -cy, cy])
drawnow
end

function F = incompleteEllipticFunction(x,k)
F = ellipticF(asin(x),k);
end

function K = completeEllipticFunction(k)
K = incompleteEllipticFunction(1,k);
end

function f = implicitMequation(x,T)
f = completeEllipticFunction(x)*T - completeEllipticFunction(1-x);
end

function x = computeX(ax,Kx,mx,t)
factor = ax/Kx;
xp = sqrt(1-t);
x = factor*incompleteEllipticFunction(xp,mx);
end

function y = computeY(ay,Ky,my,M,t)
factor = ay/Ky;
yp = sqrt(1-M/t);
y = factor*incompleteEllipticFunction(yp,my);
end
