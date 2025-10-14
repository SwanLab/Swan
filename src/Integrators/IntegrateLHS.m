function LHS = IntegrateLHS(f,test,trial,mesh,type,quadOrder)
    if nargin < 6 || isempty(quadOrder)
        qTe = test.getOrderNum();
        qTr = trial.getOrderNum();
        quadOrder = qTe + qTr;
    end
    switch type
        case 'Domain'
            lhs = integrateElementalLHS(f,test,trial,mesh,quadOrder);
            LHS = assembleMatrix(lhs,test,trial);
        case 'Boundary'
            [bMesh, l2g] = mesh.createSingleBoundaryMesh();
            [bTest,bTrial,iGlob,jGlob] = restrictTestTrialToBoundary(bMesh,test,trial,l2g);
            lhsLoc = IntegrateLHS(f,bTest,bTrial,bMesh,'Domain',quadOrder);
            [iLoc,jLoc,vals] = find(lhsLoc);
            LHS = sparse(iGlob(iLoc),jGlob(jLoc),vals, test.nDofs,trial.nDofs);
    end
end

function lhs = integrateElementalLHS(f,test,trial,mesh,quadOrder)
    quad = Quadrature.create(mesh,quadOrder);
    xV = quad.posgp;
    w  = quad.weigp;
    lhs    = zeros(test.nDofsElem,trial.nDofsElem,mesh.nelem);
    detJ   = DetJacobian(mesh);
    v = @(i) Test(test,i);
    u = @(j) Test(trial,j);
    for i = 1:test.nDofsElem
        for j = 1:trial.nDofsElem
            int = (f(u(j),v(i)).*detJ)*w';
            lhs(i,j,:) = lhs(i,j,:) + int.evaluate(xV);
        end
    end
end

function A = assembleMatrix(Aelem,f1,f2)
    dofsF1 = f1.getDofConnec();
    if isequal(f1, f2)
        dofsF2 = dofsF1;
    else
        dofsF2 = f2.getDofConnec();
    end
    nDofs1     = numel(f1.fValues);
    nDofs2     = numel(f2.fValues);
    ndofsElem1 = size(Aelem, 1);
    ndofsElem2 = size(Aelem, 2);
    
    [iElem, jElem] = meshgrid(1:ndofsElem1, 1:ndofsElem2);
    iElem = iElem(:);
    jElem = jElem(:);
    
    dofsI = dofsF1(:, iElem);
    dofsJ = dofsF2(:, jElem);
    
    rowIdx = dofsI(:);
    colIdx = dofsJ(:);
    Aval   = permute(Aelem,[3 2 1]);
    values = Aval(:);
    A = sparse(rowIdx, colIdx, values, nDofs1, nDofs2);
end

function [bTest, bTrial, iGlob, jGlob] = restrictTestTrialToBoundary(bMesh, test, trial, l2g)
    lastDofs = (l2g * test.ndimf)';
    l2g_dof = zeros(length(lastDofs),test.ndimf);
    for i = 1:test.ndimf
        l2g_dof(:,i) = lastDofs - (test.ndimf-i);
    end
    [bTest, iGlob] = restrictFunc(bMesh,test,l2g_dof);
    [bTrial,jGlob] = restrictFunc(bMesh,trial,l2g_dof);
end

function [bFunc, gFunc] = restrictFunc(bMesh,func,l2g_map)
    if func.mesh.kFace == 0
        bFunc = func.restrictBaseToBoundary(bMesh,l2g_map);
        l2g_map = reshape(l2g_map',[],1);
        gFunc = @(iLoc) l2g_map(iLoc);
    else
        bFunc = func;
        gFunc = @(iLoc) iLoc;
    end
end