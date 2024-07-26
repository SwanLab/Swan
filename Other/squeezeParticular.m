function squeezed = squeezeParticular(varargin) % t, dim, (TOTALDIM)
    t = varargin{1};
    dim = varargin{2};
    if isa(t, 'DomainFunction')
        z = varargin{3};
        index = setdiff(1:length(z),dim);
        s.operation = @(xV) reshape(t.operation(xV),[z(index) z(dim)]);
        squeezed = DomainFunction(s);
    else
        z = size(t);
        index = setdiff(1:length(z),dim);
        squeezed = reshape(t,[z(index) z(dim)]);
    end
end