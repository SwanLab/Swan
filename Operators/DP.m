function dom = DP(A,B)
s.operation = @(xV) evaluate(A,B,xV);
dom         = DomainFunction(s);
end

function aDb = evaluate(a,b,xV)
aEval = a.evaluate(xV);
bEval2 = b.evaluate(xV);
if length(size(bEval2)==3)
bEval(:,1,:,:) = bEval2;
end
aEval = pagetranspose(aEval);
aDb = pagemtimes(aEval,bEval);
aDb = squeezeParticular(aDb,1);
end