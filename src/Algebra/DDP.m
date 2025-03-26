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
    op = DP(A,B);
    fVR = op.evaluate(xV);
end

% function fVR = evaluate(A,B,xV)
%     aEval = computeLeftSideEvaluation(A,xV);
%     bEval = computeRightSideEvaluation(B,xV);
%     bTranspose = pagetranspose(bEval);
%     fVR   = trace(pagemtimes(aEval,bTranspose));
%     fVR   = squeezeParticular(fVR,1);
% end

function aEval = computeLeftSideEvaluation(A,xV)
    aEval    = A.evaluate(xV);
    isTensor = checkTensor(A,aEval);
    if ~isTensor
        error('Not enough dimensions to contract')
    end
end

function bEval = computeRightSideEvaluation(B,xV)
    bEval    = B.evaluate(xV);
    isTensor = checkTensor(B,bEval);
    if ~isTensor
        error('Not enough dimensions to contract')
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