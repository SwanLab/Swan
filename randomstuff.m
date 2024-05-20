clc; clear; close all

X = [0, 0;
    1,0;
    1,1;
    0, 1];
connec = [1 2 3 4];

s.coord = X;
s.connec = connec;
msh = Mesh.create(s);

x = [0,0;
    1.2,0;
    1.2,0.9;
    0, 0.9];

a.fValues = x - X;
a.mesh = msh;
a.order = 'P1';
p1 = LagrangianFunction(a);

quad = Quadrature.create(msh, 1);
xG = quad.posgp;
dNdx = p1.evaluateCartesianDerivatives(xG);
Grad(p1).evaluate(xG);


nPoints  = size(xG,2);
nElem = msh.nelem;
nDimG = msh.ndim;
nDimf = p1.ndimf;

GradU = reshape(Grad(p1).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
F = eye(2) + GradU;

% Bonet, p. 197 - We need to transpose F!
F = permute(F, [2 1])
C = F'.*F;

for I = 1:2
    for J = 1:2
        C(I,J) = sum(F(:,I).*F(:,J));
    end
end

jac = det(F);
invFt = inv(F');

mu = 1;
lambda = 2;

piola = mu*(F-invFt) + lambda*log(jac).*invFt;

%% 

mu = 1;
lambda = 1;
GradU2 = ActualGrad(p1);
I33 = Identity(p1);
F = GradU2 + I33;
b = F*F';

jac = Det(F);
sigma = mu.*(b-Identity(p1))./jac + lambda.*(log(jac)).*Identity(p1)./jac;