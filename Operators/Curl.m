function dom = Curl(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf = 1;
    dom = DomainFunction(s);
end

function curl = evaluate(u, xV)
    dNdx    = u.evaluateCartesianDerivatives(xV);
    nPoints = size(xV, 2);
    uV = u.getValuesByElem();
    uP = createOrthogonal(uV);
    uP = repmat(uP,[1 1 1 nPoints]);
    uP = permute(uP,[1 2 4 3]);
    grad  = pagemtimes(dNdx,uP);
    divUp = squeeze(sum(grad,1));
    curl  = divUp;
end

function uP = createOrthogonal(u)
uP = zeros(size(u));
uP(:,1,:) = u(:,2,:);
uP(:,2,:) = -u(:,1,:);
end

