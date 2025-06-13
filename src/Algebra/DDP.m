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
    dom        = DomainFunction(s);
end

function fVR = evaluate(A,B,dimA,dimB,xV)
    aEval = A.evaluate(xV);
    bEval = B.evaluate(xV);
    ndimsA = ndims(aEval)-2; %2 for nGaus and nElem
    ndimsB = ndims(bEval)-2; %2 for nGaus and nElem
    if isempty(dimA)
        dimA = [ndimsA-1 ndimsA];
    end
    if isempty(dimB)
        dimB = [ndimsB-1 ndimsB];
    end
    fVR = pagetensorprod(aEval,bEval,dimA,dimB,ndimsA,ndimsB);
    
    if ndims(fVR) <=2
        fVR = reshape(fVR,[1 size(fVR)]);
    end
end