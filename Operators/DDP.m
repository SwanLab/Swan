function dom = DDP(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    s.ndimf     = A.ndimf;
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,B,xV)
    aEval = computeLeftSideEvaluation(A,xV);    
    bEval = computeRightSideEvaluation(B,xV);
    if ndims(aEval) == ndims(bEval)
        AddB  = pagemtimes(aEval,bEval);
        fVR   = squeezeParticular(AddB, 2);
    else
        fVR = zeros(size(aEval,1),size(bEval,1),size(bEval,3),size(bEval,4));
        for idim = 1:size(aEval,1)
            a = aEval(idim,:,:,:,:);
            a = squeezeParticular(a,1);
            AddB  = pagemtimes(a,bEval);
            fVR(idim,:,:,:)   = squeezeParticular(AddB, 2);
        end
    end
end

function aEval = computeLeftSideEvaluation(A,xV)
    res      = A.evaluate(xV);
    n        = ndims(res);
    isTensor = n>=4;
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
    isTensor = n>=4;
    switch isTensor
        case true
            bEval = res;
        otherwise
            bEval(:,1,:,:) = res;
    end
end