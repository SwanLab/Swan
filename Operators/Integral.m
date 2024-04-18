function fVR = Integral(operation,funs,mesh,quad)
    xV = quad.posgp;
    dVolu(1,1,:,:)  = mesh.computeDvolume(quad);
    integ = bsxfun(@times, operation.evaluate(xV), dVolu);
    integ = squeezeParticular(sum(integ,3),3);

    fVR = assemble(integ, funs);
end

function A = assemble(Aelem, funs)
    % funs: test, trial. in that order.
    nfuns = length(funs);
    switch nfuns
        case 1
            A = assembleVector(Aelem, funs(1));

        case 2
            A = assembleMatrix(Aelem, funs(1), funs(2));

        otherwise
            error('Unknown number of functions to assemble')
    end
end


function A = assembleMatrix(Aelem, test, trial)
    dofsF1 = test.getConnec();
    if isequal(test, trial)
        dofsF2 = dofsF1;
    else
        dofsF2 = trial.getConnec();
    end
    
    nElem = size(Aelem,3);
    nDofs1 = numel(test.fValues);
    nDofs2 = numel(trial.fValues);
    ndofsElem1 = size(Aelem,1);
    ndofsElem2 = size(Aelem,2);

    res = zeros(ndofsElem1*ndofsElem2 * nElem, 3);
    strt = 1;
    fnsh = nElem;
    for i = 1:ndofsElem1
        dofsI = dofsF1(:,i);
        for j = 1:ndofsElem2
            dofsJ = dofsF2(:,j);
            a = squeeze(Aelem(i,j,:));
            matRes = [dofsI, dofsJ, a];
            res(strt:fnsh,:) = matRes;
            strt = strt + nElem;
            fnsh = fnsh + nElem;
        end
    end
    A = sparse(res(:,1), res(:,2), res(:,3), nDofs1, nDofs2);

end

function V = assembleVector(F, test)
    % Via indices
    dofConnec = test.getConnec();
    nDofsEl   = size(dofConnec,2);
    nDofs     = max(max(dofConnec));
    nGaus     = size(F,2);
    nElem     = size(F,3);
    strt = 1;
    fnsh = nElem;
    res = zeros(nDofsEl * nElem, 2);
    for iDof = 1:nDofsEl
        for igaus = 1:nGaus
            dofs = dofConnec(:,iDof);
            c = squeeze(F(iDof,igaus,:));
            matRes = [dofs, c];
            res(strt:fnsh,:) = matRes;
            strt = strt + nElem;
            fnsh = fnsh + nElem;
        end
    end
    V = sparse(res(:,1), 1, res(:,2), nDofs, 1);
end

% function dom = Integral(operation,mesh)
%     s.operation = @(quad) evaluate(operation, mesh, quad);
%     dom = IntegralFunction(s);
% end
% 
% function fVR = evaluate(operation, mesh, quad)
%     xV = quad.posgp;
%     dVolu(1,1,:,:)  = mesh.computeDvolume(quad);
%     integ = bsxfun(@times, operation.evaluate(xV), dVolu);
%     integ = squeezeParticular(sum(integ,3),3);
%     fVR = integ;
% end