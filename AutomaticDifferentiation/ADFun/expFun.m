function h = expFun(u)

if ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(exp(u), 0);
else
    h = ValDerForward(exp(u.val), exp(u.val)*u.der);
end

end