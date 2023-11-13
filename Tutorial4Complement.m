
clear;

% Mesh
x1 = linspace(0,1,10);
x2 = linspace(0,1,10);
[xv,yv] = meshgrid(x1,x2);
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh = Mesh(m);

% Previous definitions
fP1    = P1Function.create(mesh,1);
fP0    = P0Function.create(mesh,1);
dofsP1 = fP1.nDofs;
dofsP0 = fP0.nDofs;

% P0(N1(x))
IP1P0 = zeros(dofsP0,dofsP1);
for i = 1:dofsP1
    fV          = zeros(dofsP1,1);
    fV(i)       = 1;
    s.fValues   = fV;
    s.mesh      = mesh;
    fP1.fValues = fV;
    fP0         = fP1.project('P0');
    fVP0        = fP0.fValues;
    IP1P0(:,i)  = fVP0;
end

% P1(N0(x))
IP0P1 = zeros(dofsP1,dofsP0);
for i = 1:dofsP0
    fV          = zeros(dofsP0,1);
    fV(i)       = 1;
    s.fValues   = fV;
    s.mesh      = mesh;
    fP0.fValues = fV;
    fP1         = fP0.project('P1');
    fVP1        = fP1.fValues;
    IP0P1(:,i)  = fVP1;
end

% P1(N1(x))
IP1P1 = zeros(dofsP1,dofsP1);
for i = 1:dofsP1
    fV          = zeros(dofsP1,1);
    fV(i)       = 1;
    s.fValues   = fV;
    s.mesh      = mesh;
    fP1.fValues = fV;
    fP1         = fP1.project('P1');
    fVP1        = fP1.fValues;
    IP1P1(:,i)  = fVP1;
end

% Example
ss.fHandle = @(x) 1+x(1,:,:).^2+2*x(2,:,:).^2;
ss.ndimf   = 1;
ss.mesh    = mesh;
AFun       = AnalyticalFunction(ss);
fun        = AFun.project('P1');
funP0      = fun.project('P0');

x1 = funP0.fValues;

fV1 = fun.fValues;
x2 = IP1P0*fV1;