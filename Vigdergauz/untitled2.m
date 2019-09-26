function untitled2
% r = sym('r','real');
% t = sym('t','real');
% ax = sym('ax','real');
% ay = sym('ay','real');
% 
% eq(1) = r*(1-ax) -(1-ay);
% eq(2) = t - ax - ay + 1;
% 
eps = 8*1e-2;
n = 100;
theta(:,1) = linspace(0.0005,0.85,n);
%a = linspace(0+eps,1-eps,n);
for it = 1:length(theta)
    a(it,1) = computeA(theta(it));
    m(it,1) = computeM(a(it),a(it),1);
    M(it,1) = (1-1/m(it))^2;
    K(it,1) = computeK(m(it));
    r(it,1) = sqrt(theta(it)*4/pi);
    r2(it,1) = sqrt(1-M(it))*a(it)/K(it);
    
    nt = 100;
    t(:,1) = linspace(M(it),1,nt);
    xp = sqrt(1-t);
    yp = sqrt(1-M(it)./t);
    fact1(it,1) = a(it)/K(it);
    fact2(it,1) = a(it)/K(it)*sqrt(1-M(it));
    fact3(it,1) = a(it)/K(it)*ellipticF(sqrt(1-M(it)),m(it));
    fact4(it,1) = a(it)*ellipticF(sqrt(1-M(it)),m(it));
    fact5(it,1) = ellipticF(sqrt(1-M(it)),m(it));
    fact6(it,1) = ellipticF(sqrt(1-M(it)),m(it))/ellipticF(1,m(it));
    
ellipticF(1,m(it)) - ellipk(m(it))
%     for jt = 1:nt
%        Fx(jt,1) = ellipticF(xp(jt),m(it)); 
%        Fy(jt,1) = ellipticF(yp(jt),m(it));
%        x(jt,1) = fact1(it)*Fx(jt);
%        y(jt,1) = fact1(it)*Fy(jt);  
%        txi(jt,1) = Fy(jt,1)/K(it);
%     end
    
end
plot(theta,a)
plot(theta,m)
plot(m,K)
plot(theta,K)
plot(theta,M);
plot(theta,fact1,'-+')
plot(theta,fact2,'-+')
plot(theta,fact3,'-+')
plot(theta,[fact4*1.65 sqrt(4*theta/pi) sqrt(4*theta)],'-+')
plot(theta,fact5,'-+')
plot(theta,fact6,'-+')
plot(theta,[r r2])
end

function a = computeA(theta)
a = (1+theta)/2;
end

function my = computeM(ax,ay,cx)
f = @(x) cx- ax - ay*ellipk(1-x)/ellipk(x);
eps = 1e-16;
my = fzero(f,[eps,1-eps]);
end

function K = computeK(mx)
K = ellipk(mx);
end