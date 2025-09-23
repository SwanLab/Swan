function LHS = IntegrateLHSBoundary(f,test,trial,mesh,quadOrder)
if nargin < 5 || isempty(quadOrder)
    qTe = test.getOrderNum();
    qTr = trial.getOrderNum();
    quadOrder = qTe + qTr;
end
[bMesh, l2g] = mesh.createSingleBoundaryMesh();
ndof         = trial.nDofs;
LHS          = sparse(ndof,ndof);
bTest        = LagrangianFunction.create(bMesh,test.ndimf,test.order);
bTrial       = LagrangianFunction.create(bMesh,trial.ndimf,trial.order);
lhsLoc       = IntegrateLHS(f,bTest,bTrial,bMesh,quadOrder);
LHS(l2g,l2g) = lhsLoc;
end