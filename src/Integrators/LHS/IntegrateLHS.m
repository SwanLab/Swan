function LHS = IntegrateLHS(f,test,trial,mesh,type,quadOrder)
    if nargin < 6 || isempty(quadOrder)
        qTe = test.getOrderNum();
        qTr = trial.getOrderNum();
        quadOrder = qTe + qTr;
    end
    switch type
        case 'Domain'
            s.trial = trial;
            s.test  = test;
            s.mesh  = mesh;
            s.quadratureOrder = quadOrder;
            lhs = LHSIntegrator(s);
            LHS = lhs.compute(f);
        case 'Boundary'
            [bMesh, l2g] = mesh.createSingleBoundaryMesh();
            [bTest,bTrial,iGlob,jGlob] = restrictTestTrialToBoundary(bMesh,test,trial,l2g);
            lhsLoc = IntegrateLHS(f,bTest,bTrial,bMesh,'Domain',quadOrder);
            [iLoc,jLoc,vals] = find(lhsLoc);
            LHS = sparse(iGlob(iLoc),jGlob(jLoc),vals, test.nDofs,trial.nDofs);
    end
end

function [bTest, bTrial, iGlob, jGlob] = restrictTestTrialToBoundary(bMesh, test, trial, l2g)
    l2g_dof        = (l2g * test.ndimf)' - (test.ndimf-1:-1:0)';
    [bTest, iGlob] = restrictFunc(bMesh,test,l2g_dof);
    [bTrial,jGlob] = restrictFunc(bMesh,trial,l2g_dof);
end

function [bFunc, gFunc] = restrictFunc(bMesh,func,l2g_map)
    if func.mesh.kFace == 0
        bFunc = func.restrictBaseToBoundary(bMesh);
        gFunc = @(iLoc) l2g_map(iLoc);
    else
        bFunc = func;
        gFunc = @(iLoc) iLoc;
    end
end