function dom = DP(varargin)
    A = varargin{1}; B = varargin{2};
    dimA = 2; dimB = 1;
    if nargin > 2, dimA = varargin{3}; end
    if nargin > 3, dimB = varargin{4}; end

    s.operation = @(xV) evaluate(A,B,dimA,dimB,xV);
    if isa(A,'DomainFunction')
        s.mesh = A.mesh;
    else
        s.mesh = B.mesh;
    end
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,B,dimA,dimB,xV)
    aEval = computeLeftSideEvaluation(A,dimA,xV);
    bEval = computeRightSideEvaluation(B,dimB,xV);
    fVR = pagemtimes(aEval,bEval);
    if size(fVR,1) == 1
        fVR = squeezeParticular(fVR, 1);
    elseif size(fVR,2) == 1
        fVR = squeezeParticular(fVR, 2);
    end
end

function aEval = computeLeftSideEvaluation(A,dimA,xV)
    res      = A.evaluate(xV);
    isTensor = checkTensor(A,res);
    if isTensor
        aEval = res;
        if dimA == 1
            aEval = pagetranspose(aEval);
        end
    else
        aEval(1,:,:,:) = res;
    end
end

function bEval = computeRightSideEvaluation(B,dimB,xV)
    res      = B.evaluate(xV);
    isTensor = checkTensor(B,res);
    if isTensor
        bEval = res;
        if dimB == 2
            bEval = pagetranspose(bEval);
        end
    else
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