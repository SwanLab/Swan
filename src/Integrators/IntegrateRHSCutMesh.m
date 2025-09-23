function RHS = IntegrateRHSCutMesh(f,test,cutMesh,quadOrder)
if nargin < 4 || isempty(quadOrder)
    quadOrder = 2;
end
rhs = integrateElementalRHS(f,test,cutMesh,quadOrder);
RHS = assembleVector(rhs, test);
end

function rhs = integrateElementalRHS(f,test,cutMesh,quadOrder)
quad     = Quadrature.create(cutMesh.mesh,quadOrder);
xVLoc    = quad.posgp;
isoMesh  = obtainIsoparametricMesh(cutMesh);
xV       = isoMesh.evaluate(xVLoc);
globCell = cutMesh.cellContainingSubcell;
nElem    = test.mesh.nelem;
w        = repmat(quad.weigp',[1 1 cutMesh.mesh.nelem]);
rhs      = zeros(test.nDofsElem,nElem);
detJ     = DetJacobian(cutMesh.mesh);
fdetJ    = f.*detJ;
fdetJ    = fdetJ.evaluate(xVLoc);
v        = @(i) Test(test,i);
for i = 1:test.nDofsElem
    N = v(i).evaluate(xV);
    int = squeezeParticular(pagemtimes(N.*fdetJ,w),2);
    rhs(i,:) = rhs(i,:) + accumarray(globCell,int,[nElem,1],@sum,0)';
end
end

function m = obtainIsoparametricMesh(cutMesh)
coord      = cutMesh.xCoordsIso;
nDim       = size(coord,1);
nNode      = size(coord,2);
nElem      = size(coord,3);
msh.connec = reshape(1:nElem*nNode,nNode,nElem)';
msh.type   = cutMesh.mesh.type;
s.fValues  = reshape(coord,nDim,[])';
s.mesh     = msh;
s.order    = 'P1';
m          = LagrangianFunction(s);
end

function F = assembleVector(Felem, f)
dofConnec = f.getDofConnec();
nDofs     = numel(f.fValues);
rowIdx    = dofConnec(:);
Felem = Felem';
F = sparse(rowIdx, 1, Felem(:), nDofs, 1);
end