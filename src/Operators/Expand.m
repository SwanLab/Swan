function dom = Expand(varargin)
    if nargin == 1
        a = varargin{1};
        s.operation = @(xV) evaluate(a,xV);
    else
        a = varargin{1}; b = varargin{2};
        s.operation = @(xV) evaluate(a,b,xV);
    end

    if isa(a,'DomainFunction')
        s.mesh = a.mesh;
    else
        s.mesh = b.mesh;
    end
    s.ndimf = max(a.ndimf,b.ndimf);   
    dom         = DomainFunction(s);
end

function aEval = evaluate(varargin)
    if nargin == 2
        a = varargin{1}; xV = varargin{2};
        aEval      = a.evaluate(xV);
        dims = size(aEval);
        aEval = reshape(aEval,[dims(1), 1, dims(2:end)]);
    else
        a = varargin{1}; b = varargin{2}; xV = varargin{3};
        if ~isnumeric(a)
            if ~isnumeric(b)
                aEval      = a.evaluate(xV);
                bEval      = b.evaluate(xV);
                isTensorA  = checkTensor(a,aEval);
                isTensorB  = checkTensor(b,bEval);
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