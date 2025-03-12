clear;
clc;

% Preliminaries

mesh = UnitTriangleMesh(50,50);
N    = 2;
E    = 1;
nu   = 1/3;

u1      = linspace(0,1.7,mesh.nnodes);
u2      = linspace(0,2.3,mesh.nnodes);
uValues = [u1',u2'];
u       = LagrangianFunction.create(mesh,N,'P1');
u.setFValues(uValues);
strain = SymGrad(u);
str2   = strain.evaluate([0;0]);
str2   = str2(:,:,1,1);
strain = Voigt(strain);
str1   = strain.evaluate([0;0]);
str1   = str1(:,:,1);


% Voigt form

young   = ConstantFunction.create(E,mesh);
poisson = ConstantFunction.create(nu,mesh);

s.type    = 'ISOTROPIC';
s.ptype   = 'ELASTIC';
s.ndim    = mesh.ndim;
s.young   = young;
s.poisson = poisson;
tensor    = Material.create(s);

resElem = tensor.evaluate([0;0]);
C1    = resElem(:,:,1,1);

stress = DDP(tensor,strain);
sig1   = stress.evaluate([0;0]);
sig1   = sig1(:,1,1);


% Fourth order form

mu     = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E,nu);
k      = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E,nu,N);
lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(mu,k,N);

I   = eye(N);
IkI = tensorprod(I,I);
I4  = double((reshape(1:N, [N, 1, 1, 1]) == reshape(1:N, [1, 1, N, 1])) & ...
     (reshape(1:N, [1, N, 1, 1]) == reshape(1:N, [1, 1, 1, N])));

C2 = 2*mu*I4 + lambda*IkI;

sig2 = tensorprod(C2,str2,[3,4],1:2);


% Output

disp('The elastic tensor in Voigt is:');
disp(C1);
disp('The elastic tensor without Voigt is:');
disp(C2);
disp('The strain in Voigt is:');
disp(str1);
disp('The strain without Voigt is:');
disp(str2);
disp('The stress in Voigt is:');
disp(sig1);
disp('The stress without Voigt is:');
disp(sig2);

%% C*grad(u) test
E = 1; nu = 0.3; N=2;
mesh = TriangleMesh(1,1,2,2);
u    = LagrangianFunction.create(mesh,2,'P1');
young   = ConstantFunction.create(E,mesh);
poisson = ConstantFunction.create(nu,mesh);

% Voigt 
s.type    = 'ISOTROPIC';
s.ptype   = 'ELASTIC';
s.ndim    = mesh.ndim;
s.young   = young;
s.poisson = poisson;
Cvoigt    = Material.create(s);

s.type     = 'ElasticStiffnessMatrix';
s.mesh     = mesh;
s.test     = u;
s.trial    = u;
s.material = Cvoigt;
s.quadratureOrder = 2;
LHS = LHSIntegrator.create(s);
K = full(LHS.compute());

% No voigt
mu     = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E,nu);
k      = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E,nu,N);
lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(mu,k,N);
I   = eye(N);
IkI = tensorprod(I,I);
I4  = double((reshape(1:N, [N, 1, 1, 1]) == reshape(1:N, [1, 1, N, 1])) & ...
     (reshape(1:N, [1, N, 1, 1]) == reshape(1:N, [1, 1, 1, N])));
C = 2*mu*I4 + lambda*IkI;


s.quadratureOrder;
quad = Quadrature.create(mesh, s.quadratureOrder);

nnodeElem = mesh.nnodeElem;
xV = quad.posgp;
dNdx = u.evaluateCartesianDerivatives(xV);
lhs = zeros(nnodeElem*2,nnodeElem*2,mesh.nelem);
dV   = mesh.computeDvolume(quad);
dV   = permute(dV,[3, 4, 1, 2]);
for i=1:nnodeElem
    dNdxTrial = zeros(N,N,2,4,4);
    dNdxTrial(:,1,1,:,:) = dNdx(:,i,:,:);
    dNdxTrial(:,2,2,:,:) = dNdx(:,i,:,:);
    symTrial = 0.5*(dNdxTrial + pagetranspose(dNdxTrial));
    for j=1:nnodeElem
        dNdxTest = zeros(N,N,2,4,4);
        dNdxTest(:,1,1,:,:) = dNdx(:,j,:,:);
        dNdxTest(:,2,2,:,:) = dNdx(:,j,:,:);
        symTest = 0.5*(dNdxTest + pagetranspose(dNdxTest));
        sigN = pagetensorprod(C,symTrial,[3 4],[1 2],4,3);
        Kelem = pagetensorprod(symTest,sigN,[1 2],[1 2],3,3);
        Kelem = Kelem.*dV;
        A = squeezeParticular(sum(Kelem,3),3);
        I = (2*i-1):(2*i);
        J = (2*j-1):(2*j);
        lhs(I,J,:) = lhs(I,J,:) + A;
    end
end

for i=1:4
    for j=1:4
        res(1,i,j) = tensorprod(symTest(:,:,i,j),sigN(:,:,i,j),[1 2],[1 2]);
    end
end


%%
            xV   = obj.quadrature.posgp;
            dNdx = obj.test.evaluateCartesianDerivatives(xV);
            B    = obj.computeB(dNdx);            
            dV   = obj.mesh.computeDvolume(obj.quadrature);
            dV   = permute(dV,[3, 4, 2, 1]);
            Cmat = obj.material.evaluate(xV);
            Cmat = permute(Cmat,[1 2 4 3]);
            Bt   = permute(B,[2 1 3 4]);
            BtC  = pagemtimes(Bt, Cmat);
            BtCB = pagemtimes(BtC, B);
            lhs  = sum(BtCB .* dV, 4);
