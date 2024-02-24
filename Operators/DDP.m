function dom = DDP(C, eps)
    s.operation = @(xV) evaluate(C, eps, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(C, eps, xV)
    strn = eps.evaluate(xV);
    Cv   = C.evaluate(xV);

    strs = pagemtimes(Cv,strn);
    % fVR = permute(strs, [1 3 4 2]);
    fVR = strs;
end