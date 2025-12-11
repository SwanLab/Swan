function RHS = IntegrateRHS(f,test,mesh,type,quadOrder)
sample = f(Test(test,1));
switch class(sample)
    case 'UnfittedFunction'
        s.mesh     = sample.unfittedMesh;
        s.quadType = quadOrder;
        int        = RHSIntegratorUnfitted(s);
        RHS        = int.computeUnfitted(f,test);
    case 'UnfittedBoundaryFunction'
        s.mesh     = sample.unfittedMesh;
        s.quadType = quadOrder;
        int        = RHSIntegratorUnfitted(s);
        RHS        = int.computeUnfittedBoundary(f,test);
    otherwise
        if nargin < 5 || isempty(quadOrder)
            quadOrder = 2;
        end
        switch type
            case 'Domain'
                rhs = integrateElementalRHS(f,test,mesh,quadOrder);
                RHS = assembleVector(rhs, test);
            case 'Boundary'
                [bMesh, l2g]  = mesh.createSingleBoundaryMesh();
                [bTest,iGlob] = restrictTestToBoundary(test,l2g);
                rhsLoc = IntegrateRHS(f,bTest,bMesh,'Domain',quadOrder);
                [iLoc,~,vals] = find(rhsLoc);
                RHS = sparse(iGlob(iLoc),1,vals, test.nDofs,1);
        end
end
rhs = integrateElementalRHS(f,test,mesh,quadOrder);
RHS = assembleVector(rhs, test);
end

function rhs = integrateElementalRHS(f,test,mesh,quadOrder)
quad = Quadrature.create(mesh,quadOrder);
xV = quad.posgp;
w  = quad.weigp;
rhs  = zeros(test.nDofsElem,mesh.nelem);
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

function [bTest, iGlob] = restrictTestToBoundary(test, l2g)
    lastDofs = (l2g * test.ndimf)';
    l2g_dof = zeros(length(lastDofs),test.ndimf);
    for i = 1:test.ndimf
        l2g_dof(:,i) = lastDofs - (test.ndimf-i);
    end
    [bTest, iGlob] = restrictFunc(test,l2g_dof);
end

function [bFunc, gFunc] = restrictFunc(func,l2g_map)
    if func.mesh.kFace == 0
        bFunc = func.restrictToBoundary();
        l2g_map = reshape(l2g_map',[],1);
        gFunc = @(iLoc) l2g_map(iLoc);
    else
        bFunc = func;
        gFunc = @(iLoc) iLoc;
    end
end