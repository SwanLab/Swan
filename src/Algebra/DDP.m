function dom = DDP(varargin)
    A = varargin{1}; B = varargin{2};
    dimA = []; dimB = [];
    if nargin > 2, dimA = varargin{3}; end
    if nargin > 3, dimB = varargin{4}; end

    s.operation = @(xV) evaluate(A,B,dimA,dimB,xV);
    if isa(A,'DomainFunction')
        s.mesh = A.mesh;
    else
        s.mesh = B.mesh;
    end
    if isa(A,'Material')
        nA = 0; 
    else 
        nA = A.ndimf;
    end
    if isa(B,'Material')
        nB = 0;
    else 
        nB = B.ndimf;
    end

    s.ndimf = max(nA,nB);
    dom        = DomainFunction(s);
end

function fVR = evaluate(A,B,dimA,dimB,xV)
    aEval = A.evaluate(xV);
    bEval = B.evaluate(xV);

    extraDim = computeExtraDims(A,B,xV);
    ndimsA = ndims(aEval)-extraDim; %To be adapted when ndimf is vector
    ndimsB = ndims(bEval)-extraDim; %2 for nGaus and nElem
    if isempty(dimA)
        dimA = [ndimsA-1 ndimsA];
    end
    if isempty(dimB)
        dimB = [ndimsB-1 ndimsB];
    end
    fVR = pagetensorprod(aEval,bEval,dimA,dimB,ndimsA,ndimsB);
    
    if ndims(fVR) == 2
        fVR = reshape(fVR,[1 size(fVR)]);
    end
end

function extraDim = computeExtraDims(A,B,xV)
    if any(strcmp('mesh', properties(A)))
        nelem = A.mesh.nelem;
    else
        nelem = B.mesh.nelem;
    end

    extraDim = 2;
    if nelem == 1
        extraDim = extraDim - 1;
        if size(xV,2) == 1
            extraDim = extraDim -1;
        end
    end
end