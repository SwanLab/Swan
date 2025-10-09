%% Define future inputs
input.mesh = UnitQuadMesh(2,2);
input.f1.ndimf = [2 2 2 2];
input.f1.order = 'P1';
nDofs = input.mesh.nnodes*prod(input.f1.ndimf);
input.f1.values = rand(nDofs,1);

input.f2.ndimf = [2 2];
input.f2.order = 'P1';
nDofs = input.mesh.nnodes*prod(input.f2.ndimf);
input.f2.values = rand(nDofs,1);

%% Compute current solution 
f1 = LagrangianFunction.create(input.mesh,prod(input.f1.ndimf),input.f1.order);
f1.setFValues(reshape(input.f1.values,[prod(input.f1.ndimf) input.mesh.nnodes])');
f2 = LagrangianFunction.create(input.mesh,prod(input.f2.ndimf),input.f2.order);
f2.setFValues(reshape(input.f2.values,[prod(input.f2.ndimf) input.mesh.nnodes])');
x = DDP(f1,f2);
f1Val = pagetranspose(reshape(f1.evaluate([0;0]),[input.f1.ndimf,1,input.mesh.nelem]));
% f2Val = f2.evaluate([0;0]);
% xRef = pagetensorprod(f1Val,f2Val,2,1,2,1)
f2Val = pagetranspose(reshape(f2.evaluate([0;0]),[input.f2.ndimf,1,input.mesh.nelem]));
xRef = pagetensorprod(f1Val,f2Val,[3 4],[1 2],4,2)
% xRef = reshape(xRef,[1 size(xRef)]);
save('TestDDPMTensorMatrix','input','xRef')