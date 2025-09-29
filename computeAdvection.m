function LHS = computeAdvection(f,test,trial,mesh,quadOrder)
            lhs = integrateElementalLHS(f,test,trial,mesh,quadOrder);
            LHS = assembleMatrix(lhs,test,trial);
 


end


function lhs = integrateElementalLHS(f,test,trial,mesh,quadOrder)
quad = Quadrature.create(mesh,quadOrder);
% xV = quad.posgp;
% w  = quad.weigp;
% lhs    = zeros(test.nDofsElem,trial.nDofsElem,mesh.nelem);
% detJ   = DetJacobian(mesh);
% v = @(i) Test(test,i);
% u = @(j) Test(trial,j);
% for i = 1:test.nDofsElem
%     for j = 1:trial.nDofsElem
%         int = (f(u(j),v(i)).*detJ)*w';
%         lhs(i,j,:) = lhs(i,j,:) + int.evaluate(xV);
%     end
% end

xV = quad.posgp;
shapesTest  = test.computeShapeFunctions(xV);
dNdxTr      = trial.evaluateCartesianDerivatives(xV);
fG          = squeeze(f.evaluate(xV));

dVolu = mesh.computeDvolume(quad);
nGaus = quad.ngaus;
nElem = size(dVolu,2);

nNodeTest  = size(shapesTest,1);
nNodeTrial = size(dNdxTr,2);
nDofTest   = nNodeTest*test.ndimf;
nDofTrial  = nNodeTrial*trial.ndimf;


lhs = zeros(nDofTest,nDofTrial,nElem);
%             for iGaus = 1:nGaus
%                 dV(1,1,:) = dVolu(iGaus,:)';
%                 for iDof = 1:nDofETs
%                     for jDof = 1:nDofETr
%                        Ni  = shapesTest(iDof,iGaus,:);
%                        dNj = squeeze(dNdxTr(:,jDof,:,iGaus));
%                        df  = squeeze(fG(:,iGaus,:));
%                        int(1,1,:) = sum(Ni*df.*dNj,1);
%                        lhs(iDof,jDof,:) = lhs(iDof,jDof,:) + int.*dV;
%                     end
%
%                 end
%             end


for iGaus = 1 :nGaus
    for iNode = 1:nNodeTest
        for jNode = 1:nNodeTrial
            for iDim = 1:test.ndimf
                for jDim = 1:trial.ndimf
                    fdv = fG(iGaus,:).*dVolu(iGaus,:);
                    idof = test.ndimf*(iNode-1)+iDim;
                    jdof = trial.ndimf*(jNode-1)+jDim;
                    Ni = shapesTest(iNode,iGaus,:);
                    dNj = squeeze(dNdxTr(iDim,jNode,iGaus,:));
                    v = squeeze(Ni.*dNj.*(fdv'));
                    lhs(idof, jdof, :)= squeeze(lhs(idof,jdof,:)) ...
                        + v(:);
                end
            end
        end
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