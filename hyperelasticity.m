% Material parameters
mesh = UnitHexaMesh(7,7,7);
material.lambda = 0.6;
material.mu = 1;


% Creem uFun
sAF.fHandle = @(x) [0*x(1,:,:) + x(2,:,:);
                    -0.5*x(2,:,:);
                    0*x(3,:,:)];
sAF.ndimf   = 3;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

uFun = xFun.project('P1');

s.mesh= mesh;
s.material= material;
neo = NeohookeanFunctional(s);
psi = neo.compute(uFun);

quad = Quadrature.create(mesh, 2);
xG = quad.posgp;
nPoints  = quad.ngaus;
nElem = mesh.nelem;
nDimG = mesh.ndim;
nDimf = uFun.ndimf;
GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);

I33 = zeros(size(GradU));
I33(1,1,:,:) = 1;
I33(2,2,:,:) = 1;
I33(3,3,:,:) = 1;

% NOTE: We need to transpose F!
F = I33 + GradU; % deformation gradient
F = permute(F, [2 1 3 4]); % F: nDimf, nDimG, nGaus, nElem
invF = MatrixVectorizedInverter.computeInverse(F);
invFt = permute(invF, [2 1 3 4]);

jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);

% Piola
piola = material.mu*(F-invFt) + material.lambda*(jac-1).*jac.*invFt;

% grad(deltaU)
test = LagrangianFunction.create(mesh,3,'P1');
dV(1,1,:,:) = mesh.computeDvolume(quad);
dNdx = test.evaluateCartesianDerivatives(quad.posgp); % ndimG, nnodE, nG, nEl
nNodeE = size(dNdx,2);

mult = pagemtimes(piola, dNdx);
intI = mult.*dV;

s.mesh = obj.mesh;
s.material = obj.material;
test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
rhs = RHSintegrator_FirstPiola(s);
intfor = rhs.compute(obj.uFun, test);
% idea: pillar un elasticproblem, guardar u, comparar hessian amb K.

function kron_top
end

function kron_bot
end
