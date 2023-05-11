function h = mrdivideFun(u,v)

if ~isa(u,'ValDerForward') && ~isa(v,'ValDerForward') %u and v are scalars
    h = ValDerForward(u / v, 0);

elseif ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(u / v.val, (- v.der * u.val) / v.val ^ 2);

elseif ~isa(v,'ValDerForward') %v is a scalar
    h = ValDerForward(u.val / v, u.der / v);

else
    h = ValDerForward(u.val / v.val, (u.der * v.val - v.der * u.val) / v.val ^ 2);

end
end