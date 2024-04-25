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

% Products
C = kron_top(I33,I33);
C = kron_topF(I33,I33);


% s.mesh = obj.mesh;
% s.material = obj.material;
% test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
% rhs = RHSintegrator_FirstPiola(s);
% intfor = rhs.compute(obj.uFun, test);
% idea: pillar un elasticproblem, guardar u, comparar hessian amb K.

function C = OuterProductDelta(A, B)  % version 1
    C = zeros([size(A),size(B)]);
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            for k = 1:size(B,1)
                for l = 1:size(B,2)
                    C(i,j,k,l) = A(i,j)*B(k,l);
                end
            end
        end
    end
end

function C = OuterProduct(A, B)  % version 5
    C = reshape(A(:) * B(:).', [size(A), size(B)]);
end

function C= kron_topF(A,B)
    C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
    for i = 1:size(A,1)
        for k = 1:size(B,1)
            C(i,:,k,:,:,:) = pagemtimes( A(i,k,:,:), B);
        end
    end
end

function C= kron_top(A,B)
    C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            for k = 1:size(B,1)
                for l = 1:size(B,2)
                    C(i,j,k,l,:,:) = A(i,k,:,:).*B(j,l,:,:);
                end
            end
        end
    end
end

function C= kron_botF(A,B)
    C = kron_topF(A,B);
    C = pagetranspose(C);
end

function C= kron_bot(A,B)
    % isequal(pagetranspose(kron_bot(eye(3),eye(3))), kron_top(eye(3),eye(3)))
    C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            for k = 1:size(B,1)
                for l = 1:size(B,2)
                    C(i,j,k,l,:,:) = A(i,l,:,:).*B(j,k,:,:);
                end
            end
        end
    end
end