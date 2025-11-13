function RHS = IntegrateRHSCutMesh(f,test,cutMesh,quadOrder)
if nargin < 4 || isempty(quadOrder)
    quadOrder = 2;
end
rhs = integrateElementalRHS(f,test,cutMesh,quadOrder);
RHS = assembleVector(rhs, test);
end

function rhs = integrateElementalRHS(f,test,cutMesh,quadOrder)
quad     = Quadrature.create(cutMesh.mesh,quadOrder);
xV       = quad.posgp;
globCell = cutMesh.cellContainingSubcell;
nElem    = test.mesh.nelem;
w        = quad.weigp;
rhs      = zeros(test.nDofsElem,nElem);
detJ     = DetJacobian(cutMesh.mesh);
v        = @(i) Test(test,i);
for i = 1:test.nDofsElem
    int = (f(v(i)).*detJ)*w';
    int = squeezeParticular(int.evaluate(xV),2);
    rhs(i,:) = rhs(i,:) + accumarray(globCell,int,[nElem,1],@sum,0)';
end
end

function F = assembleVector(Felem, f)
dofConnec = f.dofConnec;
nDofs     = numel(f.fValues);
rowIdx    = dofConnec(:);
Felem = Felem';
F = sparse(rowIdx, 1, Felem(:), nDofs, 1);
end