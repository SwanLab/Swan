function understundingVigdergauz
x0 = [0.8,0.8];
fun = @(x) computeEquations(x);

mesh = createMesh;
% problem.objective = fun;
% problem.x0 = x0;
% problem.solver = 'fsolve';
% problem.lb = [0,0];
% problem.ub = [1,1];
% problem.options = optimset('Display','iter');
%
% [x,fsol] = fsolve(problem)

r2o = 0.998709640937773;
r1o = 0.998709640937773;

r0 = [r1o,r2o];
res1 = computeTxi(r0);
res2 = computeRho(r0);

eps = 0.05;

r1 = 0.97;
r2 = 0.99;%z2new = r2;


TOL = 1e-12;

%[r,error1] = solveCase1(r1,r2,TOL);
[r,error2] = solveCase2(r1,r2,TOL,mesh);
%[r,error3] = solveCase3(r1,r2,TOL);




figure(1)
semilogy(error1(1:2:end,1))
hold on
semilogy(error2(1:2:end,1))
%hold on
%semilogy(error3(1:2:end,1))

end

function m = createMesh()
[coord,connec] = readFile();
s.coord = coord(:,2:3);
s.connec = connec(:,2:4);
m = Mesh.create(s);
end

function [coord,connec] = readFile()
run('test2d_micro');
end




function [z,error] = solveCase1(z1,z2,TOL)
error = 1;
i = 1;
while error(i) > TOL
    [x1new,x2new] = resolvent1(z1,z2);
    
    z1new = 0.5*(x1new + z1);
    z2new = 0.5*(x2new + z2);
    
    
    i = i+1;
    error(i,1) = computeError([z1new,z2new]);
    
    [x1newnew,x2newnew] = resolvent2(z1new,z2new);
    
    z1newnew = 0.5*(x1newnew + z1new);
    z2newnew = 0.5*(x2newnew + z2new);
    
    z1 = z1newnew;
    z2 = z2newnew;
    
    i = i+1;
    error(i,1) = computeError([z1,z2]);
    
    
end
z = [z1,z2];


end

function [z,error] = solveCase2(z1,z2,TOL,m)
error = 1;
i = 1;

z = [z1,z2];
zT(1,:) = z;
plotImages(zT,z,error,m)



a = 2;
b = -1;
eps = 1e-12;
while error(i) > TOL
    
    
    xNew = resolvent2(z);
    %     z1new = z1;
    %     z2new = z2;
    
    %zNew(:) = min(1-eps,max(0.8,a*xNew + b*z));
    
%     c = 0.5;
%     cNew = 0.5;
    
   % c = 0;
   % cNew = 1;
%     
     c = -1;
     cNew = 2;
    
    zNew = (cNew*xNew + c*z);

%     z1new = x1new;
%     z2new = x2new;
    
    i = i+1;    
    zT(i,:) = zNew;    
    error(i,1) = computeError(zNew);    
    plotImages(zT,zNew,error,m)
    
    
    xNewNew = resolvent1(zNew);
    
    cNew = 0;
    cNewNew = 1;  
   
%     cNew = 0.5;
%     cNewNew = 0.5;    
%     
%     cNew = -1;
%     cNewNew = 2;
    
    zNewNew = cNewNew*xNewNew + cNew*zNew;
    %zNewNew(:) = min(1-eps,max(0.8,a*xNewNew + b*zNew));
    
    
%     z1newnew = (x1newnew );
%     z2newnew = (x2newnew );
    
    
    i = i+1;
    zT(i,:) = zNewNew;
    error(i,1) = computeError(zNewNew);    
    plotImages(zT,zNewNew,error,m)
    
%     c = 0;
%     cNew = 0;
%     cNewNew = 1;
    
    c = 1;
    cNew = 0;
    cNewNew = 1;
    d = -1;
        
    z = cNewNew*zNewNew + cNew*zNew + c*z + d*xNew;
    
end

z = [z1,z2];

end

function plotImages(zT,zNew,error,m)

figure(1)
delete(gca)
plot(zT(:,1),zT(:,2),'-+')
axis([0.8 1 0.8 1])

%error(1:2:end)

figure(2)
semilogy(error(1:end,1))
drawnow

ls = computeLevelSet(m.coord,zNew);

s.meshBackground = m;
s.unfittedType = 'INTERIOR';

figure(3)
delete(gca)
uM = UnfittedMesh(s);
uM.compute(ls);
uM.plot()
drawnow
end

function [z,error] = solveCase3(z1,z2,TOL)
error = 1;
i = 1;
z = [z1,z2]';
while error(i) > TOL
    
    xnew = projection2(z);
    znew = 0.5*(xnew + z);
    
    i = i+1;
    error(i,1) = computeError(znew);
    
    [xnewnew] = projection1(znew);
    
    znewnew = 0.5*(xnewnew + znew);
    
    z = znewnew;
    
    i = i+1;
    error(i,1) = computeError(z);
    
end
end


function xNew = resolvent2(x)
x1 = x(1);
x2 = x(2);
func2 = @(z2) computeTxi([x1,z2]);
x2new = solveF(func2,x2);
x1new = x1;
xNew = [x1new,x2new];
end

function xNew = resolvent1(x)
x1 = x(1);
x2 = x(2);
func = @(z1) computeRho([z1,x2]);
x1new = solveF(func,x1);
x2new = x2;
xNew = [x1new,x2new];
end

function xnew = projection2(x)
func = @(z2) computeTxi([x(1), z2]);
xnew(2,1) = projection(func,x(2));
xnew(1,1) = x(2);
end

function xnew = projection1(x)
func = @(z1) computeRho([z1, x(2)]);
xnew(1,1) = projection(func,x(1));
xnew(2,1) = x(1);
end

function [x,f] = projection(myFunc,x0)
% xB = linspace(0.01,1,100);  % Interval To Evaluate Over
% x = sort([xB,x0]);
% fx = myFunc(x);
% fxC = circshift(fx,-1,2);% Function Evaluated Over ‘x’
% cs = fx.*fxC;        % Product Negative At Zero-Crossings
%
%
% ind = find(cs(1:end-1) <= 0);
% xc = x(ind);   % Values Of ‘x’ Near Zero Crossings
%
% [~,it] = min(abs(fx(ind)));
%
% x0 = [x(ind(it)),x(ind(it) + 1)];

problem.objective = @(x) myFunc(x);
problem.x0 = x0;
problem.solver = 'fsolve';
TOL = 1e-12;
problem.options = optimset('Display','iter','TolFun',TOL);

[x,f] = fsolve(problem);
x = real(x);
end




function er = computeError(r)
err(1) = computeRho(r);
err(2) = computeTxi(r);
er = norm(err);
end

function [x,f] = solveF(func,x0)

%[ub,lb] = findRbounds(x0,func);
%lb = max(0,lb);
%ub = min(1,ub);
eps = 1e-13;
xmin = 0+eps;
xmax = 1-eps;

%[lb,ub] = findBounds2(func,xmax,xmin,x0);
lb = xmin;
ub = xmax;

problem.objective = @(x)abs(func(x));
%problem.x0 = x0;
problem.solver = 'fminbnd';
problem.x1 = lb;
problem.x2 = ub;
%problem.options = optimset('Display','iter','TolFun',1e-15);
TOL = 1e-12;
problem.options = optimset('Display','iter','TolFun',TOL,'TolX',1e-14);

[x,f] = fminbnd(problem);
end

function [lb,ub] = findBounds2(myFunc,xmax,xmin,x0)
xB = linspace(xmin,xmax,100);  % Interval To Evaluate Over
x = sort([xB,x0]);
fx = myFunc(x);
fxC = circshift(fx,-1,2);% Function Evaluated Over ‘x’
cs = fx.*fxC;        % Product Negative At Zero-Crossings

ind = find(cs(1:end-1) <= 0);
if isempty(ind)
    ind = true(size(cs));
end

xc = x(ind);   % Values Of ‘x’ Near Zero Crossings
[~,it] = min(abs(fx(ind)));

x0 = [x(ind(it)),x(ind(it) + 1)];
lb = x(ind(it));
ub = x(ind(it) + 1);
end


function [rub,rlb] = findRbounds(r0,F)
F0 = F(r0);
eps = 10^(-12);
if F0 >= 0
    r1 = r0 - eps;
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
    r1 = r0 - eps;
    F1 = F(r1);
    while F1 <= 0
        rnew = newPointBySecant(r0,r1,F0,F1);
        r0 = r1;
        F0 = F1;
        r1 = min(1,rnew);
        F1 = F(r1);
    end
    rub = r0;
    rlb = r1;
end
end


function x2 = newPointBySecant(x0,x1,f0,f1)
x2 = x1 - (x1-x0)/(f1 - f0)*f1;
end

function ls = computeLevelSet(coord,r)
y1 = coord(:,1) - 0.5;
y2 = coord(:,2) - 0.5;
[f1,f2,m1,m2,R] = computeParameters(r);

ye1 = ellipj(y1*f1/(m1/2),r(1));
ye2 = ellipj(y2*f2/(m2/2),r(2));

ls(:,1) = (1-ye1.^2).*(1-ye2.^2) - R;

isSmallerM1 = abs(y1) < m1/2;
isSmallerM2 = abs(y2) < m2/2;

isIn = isSmallerM1 & isSmallerM2;
isOut = ~isIn;
ls(isOut) = -abs(ls(isOut));
end

function [f1,f2,m1,m2,R] = computeParameters(r)
r1 = r(1);
r2 = r(2);
R  = (1-r1)*(1-r2)/(r1*r2);

F = @(x,k) ellipticF(asin(x),k);
f = @(r) F(sqrt(1-R),r);

s = @(r) f(r)/F(1,r);
t = @(r) F(1,1-r)/F(1,r);

f1 = f(r1);
f2 = f(r2);

s1 = s(r1);
s2 = s(r2);

t1 = t(r1);
t2 = t(r2);

m1 = (1-t1)/(1-t1*t2)*s1;
m2 = (1-t2)/(1-t1*t2)*s2;


end

function [res,dres] = computeTxi(r)
r = real(r);

r1 = r(1);
r2 = r(2:end);
tanXiV = tan(20*pi/80);
%tanXiV = 1;
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
res = (tanXiV - tanXi);
dres = [];
end

function [res,dres] = computeRho(r)
r = real(r);
rhoV = 0.2;
r1 = r(1:end-1);
r2 = r(end);

F = @(x,k) ellipticF(asin(x),k);
t = @(r) F(1,1-r)./F(1,r);

t1 = t(r1);
t2 = t(r2);


rho = 1 - (t1.*(1-t2) + t2.*(1 - t1))./(1-t1.*t2);
res = (rho - rhoV);
%res = real(res);
dres = [];
end
