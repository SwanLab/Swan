function A = Expand(a)
    if isa(a,'BaseFunction')
        s.operation = @(xV) evaluate(a,xV);
        s.mesh = a.mesh;
        s.ndimf = a.ndimf;
        A = DomainFunction(s);
    else
        A = a;
    end


end

function aEval = evaluate(a,xV)
    aEval      = a.evaluate(xV);
    isTensorA  = checkTensor(a,aEval);
    if ~isTensorA
        dims = size(aEval);
        aEval = reshape(aEval,[dims(1), 1, dims(2:end)]);
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