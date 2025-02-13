function dom = Expand(A)
    s.operation = @(xV) evaluate(A,xV);
    if isa(A,'DomainFunction')
        s.mesh = A.mesh;
    else
        s.mesh = B.mesh;
    end
    dom         = DomainFunction(s);
end

function aEval = evaluate(a,xV)
    aEval = adaptSize(a,xV);
end

function aEval = adaptSize(a,xV)
    aEval = a.evaluate(xV);
    if ndims(aEval)<=4
        dims = size(aEval);
        aEval = reshape(aEval,[dims(1), 1, dims(2:end)]);
    end
end