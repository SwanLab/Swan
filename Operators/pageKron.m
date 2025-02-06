function dom = pageKron(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    s.ndimf     = A.ndimf;
    s.mesh      = A.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,B,xV)
    aEval = computeLeftSideEvaluation(A,xV);
    bEval = computeRightSideEvaluation(B,xV);
    sizeA = size(aEval);
    sizeB = size(bEval);
    if sizeA(end:end) == sizeB(end:end)
        res = zeros([sizeA(1)*sizeB(1),sizeA(2)*sizeB(2),sizeA(3:end)]);
        for i=1:sizeA(1)
            for j=1:sizeA(2)
                for k=1:sizeB(1)
                    for l=1:sizeB(2)
                        indexA = [i,j,repmat({':'}, 1, ndims(aEval)-2)];
                        Aij = squeeze(aEval(indexA{:}));
                        indexB = [k,l,repmat({':'}, 1, ndims(bEval)-2)];
                        Bkl = squeeze(bEval(indexB{:}));

                        idx1 = sizeB(1)*(i-1) + k;
                        idx2 = sizeB(2)*(j-1) + l;
                        indexRes = [idx1,idx2,repmat({':'}, 1, ndims(res)-)];
                        res(indexRes{:}) = Aij.*Bkl;
                    end
                end
            end
        end
        fVR = res;
    else
        error('Array dimensions are not the same')
    end
end

function aEval = computeLeftSideEvaluation(A,xV)
    res      = A.evaluate(xV);
    n        = ndims(res);
    isTensor = n>=3;
    switch isTensor
        case true
            aEval = res;
        otherwise
            aEval(1,:,:,:) = res;
    end
end

function bEval = computeRightSideEvaluation(B,xV)
    res      = B.evaluate(xV);
    n        = ndims(res);
    isTensor = n>=3;
    switch isTensor
        case true
            bEval = res;
        otherwise
            bEval(:,1,:,:) = res;
    end
end