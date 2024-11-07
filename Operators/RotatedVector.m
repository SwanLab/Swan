function dom = RotatedVector(A,u)
    s.operation = @(xV) evaluate(A, u, xV);
    dom = DomainFunction(s);
end

function fAv = evaluate(A, u, xV)

uEval = u.evaluate(xV);
fAv   = pagemtimes(A,uEval);

end