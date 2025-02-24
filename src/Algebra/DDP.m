function dom = DDP(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    if isa(A,'DomainFunction')
        s.mesh = A.mesh;
    else
        s.mesh = B.mesh;
    end
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,B,xV)
    aEval = computeLeftSideEvaluation(A,xV);
    bEval = computeRightSideEvaluation(B,xV);
    fVR   = pagemtimes(aEval,bEval);
    if size(fVR,1) == 1
        fVR = squeezeParticular(fVR, 1);
    elseif size(fVR,2) == 1
        fVR = squeezeParticular(fVR, 2);
    end
end

function aEval = computeLeftSideEvaluation(A,xV)
    res      = A.evaluate(xV);
    isTensor = checkTensor(A,res);
    switch isTensor
        case true
            aEval = res;
        otherwise
            aEval(1,:,:,:) = res;
    end
end

function bEval = computeRightSideEvaluation(B,xV)
    res      = B.evaluate(xV);
    isTensor = checkTensor(B,res);
    switch isTensor
        case true
            bEval = res;
        otherwise
            bEval(:,1,:,:) = res;
    end
end

function isTensor = checkTensor(A,res)
    n = ndims(res);
    if isa(A,'Material')
        isTensor = true;
    else
        if A.mesh.nelem == 1
            isTensor = n>=3;
        else
            isTensor = n>=4;
        end
    end
end