function fVR = Integral(operation,funs,mesh,quad)
    xV = quad.posgp;
    dVolu(1,1,:,:)  = mesh.computeDvolume(quad);
    integ = bsxfun(@times, operation.evaluate(xV), dVolu);
    integ = squeezeParticular(sum(integ,3),3);

    fVR = assemble(integ, funs);
end

function A = assemble(Aelem, funs)
    % funs: test, trial. in that order.
    test = funs(1);
    dofsF1 = test.getConnec();
    if length(funs) == 2
        trial = funs(2);
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