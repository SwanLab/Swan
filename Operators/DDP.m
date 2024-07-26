function dom = DDP(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    s.ndimf     = max(abs(A.ndimf - B.ndimf),1);
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,B,xV)
    aEval = computeLeftSideEvaluation(A,xV);
    bEval = computeRightSideEvaluation(B,xV);
    AddB  = pagemtimes(aEval,bEval);
    fVR   = squeezeParticular(AddB, 2);
end

function aEval = computeLeftSideEvaluation(A,xV)
    res      = A.evaluate(xV);
    is4OrderTensor = A.ndimf == 6;
    switch is4OrderTensor
        case true
            aEval = res;
        otherwise
            aEval(1,:,:,:) = res;
    end

end

function bEval = computeRightSideEvaluation(B,xV)
    res      = B.evaluate(xV);
    n        = ndims(res);
    is4OrderTensor = B.ndimf == 6;
    switch is4OrderTensor
        case true
            bEval = res;
        otherwise
            bEval(:,1,:,:) = res;
    end

end