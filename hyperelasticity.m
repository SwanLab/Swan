% Material parameters
mesh = UnitHexaMesh(7,7,7);
material.lambda = 3/4;
material.mu = 3/8;


% Creem uFun
sAF.fHandle = @(x) [x(1,:,:) + x(2,:,:);
                    -2*x(2,:,:);
                    x(3,:,:)];
sAF.ndimf   = 3;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

uFun = xFun.project('P1');

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

F = I33 + GradU; % deformation gradient

% WE NEED TO TRANSPOSE F!!


s.mesh = obj.mesh;
s.material = obj.material;
test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
rhs = RHSintegrator_FirstPiola(s);
intfor = rhs.compute(obj.uFun, test);
% idea: pillar un elasticproblem, guardar u, comparar hessian amb K.