clc; clear; close all

X = [0, 0;
    4,0;
    0,3];
connec = [1 2 3];

s.coord = X;
s.connec = connec;
msh = Mesh.create(s);

x = [2,3;
    10,3;
    10,9;];

a.fValues = x - X;
a.mesh = msh;
a.order = 'P1';
p1 = LagrangianFunction(a);

quad = Quadrature.create(msh, 1);
xG = quad.posgp;
p1.evaluateCartesianDerivatives(xG);
Grad(p1).evaluate(xG);


nPoints  = size(xG,2);
nElem = msh.nelem;
nDimG = msh.ndim;
nDimf = p1.ndimf;

GradU = reshape(Grad(p1).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
F = eye(2) + GradU

% Bonet, p. 197 - We need to transpose F!
F = permute(F, [2 1])
C = F'.*F

for I = 1:2
    for J = 1:2
        C(I,J) = sum(F(:,I).*F(:,J));
    end
end