
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



% Open Product test (gerard things)
C1op = squeeze(tensorprod(sig1,sig1));
pagemtimes(sig1,sig1')
C2op = squeeze(tensorprod(sig2,sig2));

sig1op = C1op*str1
sig2op = tensorprod(C2op,str2,[3,4],1:2)
