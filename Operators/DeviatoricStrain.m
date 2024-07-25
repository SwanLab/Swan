function dom = DeviatoricStrain(u)
    s.operation = @(xV) evaluate(u,xV);
    dom         = DomainFunction(s);
end

function fEval = evaluate(u, xV)
    e     = AntiVoigt(SymGrad(u));
    ed    = Voigt(Deviatoric(e));
    fEval = ed.evaluate(xV);

end