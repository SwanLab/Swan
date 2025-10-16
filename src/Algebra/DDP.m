function dom = DDP(varargin)
    A = varargin{1}; B = varargin{2};
    ndimsA = length(A.ndimf); ndimsB = length(B.ndimf); 
    dimA = [ndimsA-1 ndimsA]; dimB = [ndimsB-1 ndimsB];
    if nargin > 2, dimA = varargin{3}; end
    if nargin > 3, dimB = varargin{4}; end

    ndimfRes = [A.ndimf(setdiff(1:end,dimA)), B.ndimf(setdiff(1:end,dimB))];
    if isempty(ndimfRes) s.ndimf = 1; else s.ndimf = ndimfRes; end
    s.operation = @(xV) evaluate(A,B,dimA,dimB,ndimsA,ndimsB,ndimfRes,xV);
    s.mesh  = A.mesh;
    dom     = DomainFunction(s);
end

function fVR = evaluate(A,B,dimA,dimB,ndimsA,ndimsB,ndimfRes,xV)
    aEval = A.evaluate(xV);
    bEval = B.evaluate(xV);
    fVR = pagetensorprod(aEval,bEval,dimA,dimB,ndimsA,ndimsB);
    if isempty(ndimfRes)
        fVR = reshape(fVR,[1 size(fVR)]);
    end
end