function dom = RotatedVector(A,u)
    s.operation = @(xV) evaluate(A, u, xV);
    dom = DomainFunction(s);
end

function fAv = evaluate(A, u, xV)

uEval     = u.evaluate(xV);
dim       = size(uEval,1);
ngauss    = size(uEval,2);
nElem     = size(uEval,3);
uReshaped = zeros(dim,1,ngauss,nElem);

uReshaped(:,1,:,:) = uEval;

fAv = pagemtimes(A,uReshaped);
fAv = squeezeParticular(fAv,2);

end