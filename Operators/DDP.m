function dom = DDP(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,B,xV)
    aEval = computeLeftSideEvaluation(A,xV);
    bEval = computeRightSideEvaluation(B,xV);
    AddB  = pagemtimes(aEval,bEval);
    if size(AddB,1) == 1
        fVR   = squeezeParticular(AddB, 1);
    elseif size(AddB,2) == 1
        fVR   = squeezeParticular(AddB, 2);
    end
end

function aEval = computeLeftSideEvaluation(A,xV)
    res      = squeezeParticular(A.evaluate(xV),[1 2]);
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
    res      = squeezeParticular(B.evaluate(xV),[1 2]);
    n        = ndims(res);
    isTensor = n>=4;
    switch isTensor
        case true
            bEval = res;
        otherwise
            bEval(:,1,:,:) = res;
    end
end