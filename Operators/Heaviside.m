function dom = Heaviside(x)
    s.operation = @(xV) evaluate(x, xV);
    s.ndimf = x.ndimf;
    dom = DomainFunction(s);
end

function y = evaluate(x, xV)
    xEval = x.evaluate(xV);
    y = zeros(size(xEval));
    y(xEval > 0) = 1;
end