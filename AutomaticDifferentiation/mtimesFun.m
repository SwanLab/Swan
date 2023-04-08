function h = mtimesFun(u,v) %VALDER/MTIMES overloads * with at least one valder

if ~isa(u,'ValDerForward') && ~isa(v,'ValDerForward') %u and v are scalars
    h = ValDerForward(u * v, 0);

elseif ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(u * v.val, u * v.der);

elseif ~isa(v,'ValDerForward') %v is a scalar
    h = ValDerForward(v * u.val, v * u.der);

else
    h = ValDerForward(u.val * v.val, u.der * v.val + u.val * v.der);

end
end