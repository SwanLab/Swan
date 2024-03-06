function dom = DDP(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    s.ndimf     = A.ndimf/B.ndimf;
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,B,xV)
    aEval = computeLeftSideEvaluation(A,xV);
    bEval = computeRightSideEvaluation(B,xV);
    fVR = pagemtimes(aEval,bEval);
end

function aEval = computeLeftSideEvaluation(A,xV)
    res      = A.evaluate(xV);
    if size(res,2) == 1
        aEval = pagetranspose(res);
    else
        aEval = res;
    end
end

function bEval = computeRightSideEvaluation(B,xV)
    res      = B.evaluate(xV);
    if size(res,1) == 1
        bEval = pagetranspose(res);
    else
        bEval = res;
    end
end