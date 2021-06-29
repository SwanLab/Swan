function PlaneOfResiduals

rhoOptimal = 0.2;%0.61;
tanXiV = tan(20*pi/80);
%tanXiV = 0.9;


r1_0 = 0.97;
r2_0 = 0.99;

r2o = 0.998709640937773;
r1o = 0.998709640937773;

r0 = [r2o;r1o]

computeResidual(r0,rhoOptimal,tanXiV)


eps = 1e-13;
n1 = 20;
n2 = 20;
r1V = linspace(0.8,1-eps,n1);
r2V = linspace(0.8,1-eps,n2);

for i = 1:n1
    for j = 1:n2
        r1 = r1V(i);
        r2 = r2V(j);
        res = computeResidual([r1,r2],rhoOptimal,tanXiV);
        ResTxi(i,j) = abs(res(1));
        ResRho(i,j) = abs(res(2));

    end  
    i
end

figure(1)
surf(r1V,r2V,ResTxi')
figure(2)
contour(r1V,r2V,ResTxi,50)

figure(3)
surf(r1V,r2V,ResRho')
figure(4)
contour(r1V,r2V,ResRho,50)








end

function res = computeResidual(r,rhoOptimal,tanXiV)

res(1) = computeTxi(r,tanXiV);
res(2) = computeRho(r,rhoOptimal);

end


function [res,dres] = computeTxi(r,tanXiV)
r = real(r);

r1 = r(1);
r2 = r(2:end);
R  = (1-r1).*(1-r2)./(r1.*r2);

F = @(x,k) ellipticF(asin(x),k);
f = @(r) F(sqrt(1-R),r);
s = @(r) f(r)./F(1,r);
t = @(r) F(1,1-r)./F(1,r);

s1 = s(r1);
s2 = s(r2);
t1 = t(r1);
t2 = t(r2);

m1 = (1-t1)./(1-t1.*t2).*s1;
m2 = (1-t2)./(1-t1.*t2).*s2;

tanXi = (1-t1)./(1-t2).*s1./s2;
%res = tan(txiV) - tanXi;
res = tanXiV - tanXi;
dres = [];
end

function [res,dres] = computeRho(r,rhoV)
r = real(r);
r1 = r(1:end-1);
r2 = r(end);

F = @(x,k) ellipticF(asin(x),k);
t = @(r) F(1,1-r)./F(1,r);

t1 = t(r1);
t2 = t(r2);


rho = 1 - (t1.*(1-t2) + t2.*(1 - t1))./(1-t1.*t2);
res = rho - rhoV;
res = real(res);
dres = [];
end