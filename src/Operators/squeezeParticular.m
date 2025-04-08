function squeezed = squeezeParticular(input,dim) % array, dim
    if isa(input, 'BaseFunction')
        s.operation = @(xV) evaluate(input,dim,xV);
        s.mesh   = input.mesh;
        s.ndimf  = input.ndimf;
        squeezed = DomainFunction(s);
    else
        squeezed = squeeze(input,dim);
    end
end

function squeezed = evaluate(t, dim, xV)
    array = t.evaluate(xV);
    squeezed = squeeze(array,dim);
end

function squeezed = squeeze(array,dim)
    sz = size(array);
    idxNoUnit = sz~=1;
    idx2remove = false(1,length(sz));
    idx2remove(dim) = true;
    idx2keep = idxNoUnit | ~idx2remove;
    squeezed = reshape(array,[sz(idx2keep),1,1]);
end