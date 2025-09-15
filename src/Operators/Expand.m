function A = Expand(varargin)
    a = varargin{1};
    ndim = 2;
    if nargin==2
        ndim=varargin{2};
    end
    if isa(a,'BaseFunction')
        s.operation = @(xV) evaluate(a,ndim,xV);
        s.mesh = a.mesh;
        s.ndimf = a.ndimf;
        A = DomainFunction(s);
    else
        A = a;
    end

end

function aEval = evaluate(a,ndim,xV)
    aEval         = a.evaluate(xV);
    extraDims     = computeExtraDims(a,xV);
    expandTensor  = checkTensorSize(a,ndim,extraDims,aEval);
    if expandTensor
        dims = size(aEval);
        extraDimTensor = ones(1,(ndim+extraDims)-ndims(aEval));
        aEval = reshape(aEval,[dims(1), extraDimTensor, dims(2:end)]);
    end
end

function extraDim = computeExtraDims(a,xV)
    extraDim = 2;
    if a.mesh.nelem == 1
        extraDim = extraDim - 1;
        if size(xV,2) == 1
            extraDim = extraDim -1;
        end
    end
end

function expandTensor = checkTensorSize(a,ndim,extraDims,res)
    dimTensor = ndims(res);
    if isa(a,'Material')
        expandTensor = false;
    elseif dimTensor-extraDims >= ndim
        expandTensor = false;
    else
        expandTensor = true;
    end
end
