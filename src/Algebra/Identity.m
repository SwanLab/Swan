function dom = Identity(u)
    s.operation = @(xV) evaluate(u, xV);
    s.mesh = u.mesh;
    s.ndimf = u.ndimf;
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    nDimf = u.ndimf;
    nGaus = size(xV,2);
    nElem = u.mesh.nelem;
    sz = [nDimf, nDimf, nGaus, nElem];
    I = zeros(sz);
    for i = 1:nDimf
        I(i,i,:,:) = 1;
    end
    fVR = I;
end