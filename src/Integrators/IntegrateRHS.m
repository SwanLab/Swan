function RHS = IntegrateRHS(f,test,mesh,quadOrder)
if nargin < 4 || isempty(quadOrder)
    quadOrder = 2;
end
rhs = integrateElementalRHS(f,test,mesh,quadOrder);
RHS = assembleVector(rhs, test);
end

function rhs = integrateElementalRHS(f,test,mesh,quadOrder)
quad = Quadrature.create(mesh,quadOrder);
xV = quad.posgp;
w  = quad.weigp;
rhs  = zeros(test.nDofsElem,test.mesh.nelem);
detJ = DetJacobian(mesh);
v = @(i) Test(test,i);
for i = 1:test.nDofsElem
    int = (f(v(i)).*detJ)*w';
    rhs(i,:) = rhs(i,:) + squeezeParticular(int.evaluate(xV),2);
end
end

function F = assembleVector(Felem, f)
dofConnec = f.getDofConnec();
nDofs     = numel(f.fValues);
rowIdx    = dofConnec(:);
Felem = Felem';
F = sparse(rowIdx, 1, Felem(:), nDofs, 1);
end
