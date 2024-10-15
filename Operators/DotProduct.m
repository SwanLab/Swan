function dom = DotProduct(u,v)
    s.operation = @(xV) evaluate(u, v, xV);
    dom = DomainFunction(s);
end

function fAv = evaluate(u, v, xV)


ures = reshapeVector(u,xV);
vres = reshapeVector(v,xV);
ures = pagetranspose(ures);
fAv = pagemtimes(ures,vres);
fAv = squeezeParticular(fAv,2);

end

function ures = reshapeVector(u,xV)
uEval     = u.evaluate(xV);
dim       = size(uEval,1);
ngauss    = size(uEval,2);
nElem     = size(uEval,3);
ures = zeros(dim,1,ngauss,nElem);
ures(:,1,:,:) = uEval;

end