function squeezed = squeezeParticular(varargin) % t, dim, (TOTALDIM)
    t = varargin{1};
    dim = varargin{2};
    if isa(t, 'DomainFunction')
        sz = varargin{3};
        idxNoUnit = sz~=1;
        idx2remove = false(1,length(sz));
        idx2remove(dim) = true;
        idx2keep = idxNoUnit | ~idx2remove;
        s.operation = @(xV) reshape(t.operation(xV),[sz(idx2keep),1,1]);
        squeezed = DomainFunction(s);
    else
        sz = size(t);
        idxNoUnit = sz~=1;
        idx2remove = false(1,length(sz));
        idx2remove(dim) = true;
        idx2keep = idxNoUnit | ~idx2remove;
        squeezed = reshape(t,[sz(idx2keep),1,1]);
    end
end