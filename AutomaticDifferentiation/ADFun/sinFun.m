function h = sinFun(u)

if ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(sin(u), 0);
else
    h = ValDerForward(sin(u.val), cos(u.val)*u.der);
end

end