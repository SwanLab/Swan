function h = uminusFun(u)

if ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(u - v.val, v.der);

else
    h = ValDerForward(u.val - v.val, u.der - v.der);

end
end