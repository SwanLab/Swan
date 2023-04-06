function h = plusFun(u,v)

if ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(u + v.val, v.der);

elseif ~isa(v,'ValDerForward') %v is a scalar
    h = ValDerForward(v + u.val, u.der);

else
    h = ValDerForward(u.val + v.val, u.der + v.der);

end
end