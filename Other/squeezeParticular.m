function squeezed = squeezeParticular(varargin) % t, dim, (TOTALDIM)
    t = varargin{1};
    dim = varargin{2};
    if isa(t, 'DomainFunction')
        z = varargin{3};
        index = setdiff(1:length(z),dim);
        s.operation = @(xV) reshape(t.operation(xV),[z(index) z(dim)]);
        squeezed = DomainFunction(s);
    else
        sz = size(t);
        idxNoUnit = sz~=1;
        idx2remove = false(1,length(sz));
        idx2remove(dim) = true;
        idx2keep = idxNoUnit | ~idx2remove;
        squeezed = reshape(t,sz(idx2keep));
    end
end