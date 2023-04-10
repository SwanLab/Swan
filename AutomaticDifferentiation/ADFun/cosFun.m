function h = cosFun(u)

if ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(cos(u), 0);
else
    h = ValDerForward(cos(u.val), -sin(u.val)*u.der);
end

end