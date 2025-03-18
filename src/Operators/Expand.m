function dom = Expand(a,b)
    s.operation = @(xV) evaluate(a,b,xV);
    if isa(a,'DomainFunction')
        s.mesh = a.mesh;
    else
        s.mesh = b.mesh;
    end
    s.ndimf = max(a.ndimf,b.ndimf);   
    dom         = DomainFunction(s);
end

function aEval = evaluate(a,b,xV)
    if ~isnumeric(a) 
        if ~isnumeric(b)
            aEval      = a.evaluate(xV);
            bEval      = b.evaluate(xV);
            isTensorA = checkTensor(a,aEval);
            isTensorB = checkTensor(b,bEval);
            if ~isTensorA
                if isTensorB
                dims = size(aEval);
                aEval = reshape(aEval,[dims(1), 1, dims(2:end)]);
                end
            end
        else
            aEval = a.evaluate(xV);
        end
    else
        aEval = a;
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