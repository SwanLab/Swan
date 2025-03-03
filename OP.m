function dom = OP(A,B)
    s.operation = @(xV) evaluate(A,B,xV);
    if isa(A,'DomainFunction')
        s.mesh = A.mesh;
    else
        s.mesh = B.mesh;
    end
    dom         = DomainFunction(s);
end

function aDb = evaluate(a,b,xV)
    aEval = adaptSize(a,xV);
    bEval = adaptSize(b,xV);
    bEval = pagetranspose(bEval);
    aDb = pagemtimes(aEval,bEval);
    aDb = squeezeParticular(aDb,1);
end

function aEval = adaptSize(a,xV)
    aEval = a.evaluate(xV);
    if ndims(aEval)<=4
        dims = size(aEval);
        aEval = reshape(aEval,[dims(1), 1, dims(2:end)]);
    end
end