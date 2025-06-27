function dom = Shape(u)
    s.operation = @(xV) evaluate(u, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    shapes = u.computeShapeFunctions(xV);
    % shapes = repmat(shapes, [1 1 u.mesh.nelem]);
    fVR = shapes;
    % fVR = reshape(grad, [nDimG*nDimf,nPoints, nElem]);
end