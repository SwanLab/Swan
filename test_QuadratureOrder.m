o1.ngaus = 1;
o1.posgp(:,1) = [0,0];
o1.weigp = 4;

o2.ngaus = 4;
a =  0.577350269189626;
o2.posgp(:,1) = [-a,-a];
o2.posgp(:,2) = [+a,-a];
o2.posgp(:,3) = [-a,+a];
o2.posgp(:,4) = [+a,+a];
o2.weigp =  [1,1,1,1];

posgl(1) =-0.774596669241483;
posgl(2) = 0.0;
posgl(3) = 0.774596669241483;
weigl(1) = 0.555555555555556;
weigl(2) = 0.888888888888889;
weigl(3) = 0.555555555555556;
o3.ngaus = 9;
igaus = 0;
nlocs = 3;
for ilocs = 1:nlocs
    for jlocs = 1:nlocs
        igaus = igaus+1;
        o3.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
        o3.posgp(1,igaus) = posgl(ilocs);
        o3.posgp(2,igaus) = posgl(jlocs);
    end
end

posgl(1) = sqrt(3/7-2/7*sqrt(6/5));
posgl(2) =-sqrt(3/7-2/7*sqrt(6/5));
posgl(3) = sqrt(3/7+2/7*sqrt(6/5));
posgl(4) =-sqrt(3/7+2/7*sqrt(6/5));
weigl(1) = (18+sqrt(30))/36;
weigl(2) = (18+sqrt(30))/36;
weigl(3) = (18-sqrt(30))/36;
weigl(4) = (18-sqrt(30))/36;
o4.ngaus = 16;
igaus = 0;
for ilocs = 1:4
    for jlocs = 1:4
        igaus = igaus+1;
        o4.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
        o4.posgp(1,igaus) = posgl(ilocs);
        o4.posgp(2,igaus) = posgl(jlocs);
    end
end

%%
o1.ngaus = 1;
o1.weigp = [1/2];
o1.posgp = [1/3;1/3];

o2.ngaus = 3;
o2.weigp = [1/6;1/6;1/6];
o2.posgp = [2/3,1/6;1/6,1/6;1/6,2/3]';

o3.ngaus = 4;
o3.weigp = [-27/96;25/96;25/96;25/96];
o3.posgp = [1/3,1/3;1/5,1/5;3/5,1/5;1/5,3/5]';

o4.ngaus = 6;
o4.weigp = [0.223381589678011/2;0.223381589678011/2;0.223381589678011/2;0.109951743655322/2;0.109951743655322/2;0.109951743655322/2];
o4.posgp = [0.445948490915965,0.445948490915965,0.108103018168070,0.091576213509771,0.091576213509771,0.816847572980459;
    0.108103018168070,0.445948490915965,0.445948490915965,0.816847572980459,0.091576213509771,0.091576213509771];

%%

o1ngaus = 1; % tetrahedra
o1.weigp = 1/6;
o1.posgp = [1/4;1/4;1/4];

o2.ngaus = 4;
a = 0.58541020;
b = 0.13819660;
o2.posgp = [a,b,b,b;
    b,a,b,b;
    b,b,a,b];
o2.weigp = [0.041666667,0.041666667,0.041666667,0.041666667];

o3.ngaus = 5;
a = 0.25;
b = 0.5;
c = 1/6;
o3.posgp = [a,b,c,c,c;
    a,c,c,c,b;
    a,c,c,b,c];
o3.weigp = 1/6*[-0.8,0.45,0.45,0.45,0.45];

%% 
syms x y z
f(x,y) = x^3;

intAn = int(int(int(f,z,0,1-x-y),x,0,1-y),y,0,1);

oArray = [o1 o2 o3];
intNum = zeros(size(oArray));
for o = 1:length(oArray)
    for i = 1:oArray(o).ngaus
        intNum(o) = intNum(o) + oArray(o).weigp(i)*subs(subs(f,'x',oArray(o).posgp(1,i)),'y',oArray(o).posgp(2,i));
    end
end

display(intAn)
display(intNum)

%%
s.coord = [0 0 0;1 0 0;0 1 0;0 0 1];
s.connec = [1 2 3 4];
mesh = Mesh(s);

% mesh = UnitQuadMesh(3,3);

sAF.fHandle = @(x) x(1,:,:).*x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1fun = xFun.project('P1');
p2fun = xFun.project('P2');
p3fun = xFun.project('P3');

%%

n = 12;
a = -1;
b = 1;

[x,w] = legzo(n,a,b);

fprintf('obj.ngaus = %i;\n',length(w))
fprintf('obj.weigp = [')
for i = 1:length(w)
    if i == length(w)
        fprintf('%.14d',w(i))
    else
        fprintf('%.14d, ',w(i))
    end
end
fprintf('];\n')

fprintf('obj.posgp = [')
for i = 1:length(w)
    if i == length(w)
        fprintf('%.14d',x(i))
    else
        fprintf('%.14d, ',x(i))
    end
end
fprintf('];\n')


%%

n = ["ORDER4"];
for i = 1:length(n)
q{i} = Quadrature_Tetrahedra();
q{i}.computeQuadrature(n(i))
end

syms x y z
f(x,y) = (x+1)^2*(y+1)^2;

intAn = int(int(int(f,x,0,1-y-z),y,0,1-z),z,0,1);

oArray = q;
intNum = zeros(size(oArray));
for o = 1:length(oArray)
    for i = 1:oArray{o}.ngaus
        intNum(o) = intNum(o) + oArray{o}.weigp(i)*subs(subs(f,'x',oArray{o}.posgp(1,i)),'y',oArray{o}.posgp(2,i));
    end
end

display(intAn)
display(intNum)

%%

function [x,w]=legzo(n, a, b)
%       =========================================================
%       Purpose : Compute the zeros of Legendre polynomial Pn(x)
%                 in the interval [a,b], and the corresponding
%                 weighting coefficients for Gauss-Legendre
%                 integration
%       Input :   n    --- Order of the Legendre polynomial
%                 a    --- Lower boundary (optional)
%                 b    --- Upper boundary (optional)
%       Output:   x(n) --- Zeros of the Legendre polynomial
%                 w(n) --- Corresponding weighting coefficients
%       =========================================================
if nargin == 1
    a = -1;
    b =  1;
end;
x = zeros(1, n);
w = zeros(1, n);
m = (n+1)/2;
h = b-a;

for ii=1:m
    z = cos(pi*(ii-.25)/(n+.5)); % Initial estimate.
    z1 = z+1;
    while abs(z-z1)>eps
        p1 = 1;
        p2 = 0;
        for jj = 1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
        end
        pp = n*(z*p1-p2)/(z^2-1); % The L.P. derivative.
        z1 = z;
        z = z1-p1/pp;
    end
    x(ii) = z; % Build up the abscissas.
    x(n+1-ii) = -z;
    w(ii) = h/((1-z^2)*(pp^2)); % Build up the weights.
    w(n+1-ii) = w(ii);
end

if a ~= -1 || b ~= 1
    x = (x+1)*(h/2) + a;
end
end

