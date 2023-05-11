function h = mrdivideFun(u,v)

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u / v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(u / v.val, (- v.grad * u.val) / v.val ^ 2);

elseif ~isa(v,'ValGradForward') %v is a scalar
    h = ValGradForward(u.val / v, u.grad / v);

else
    h = ValGradForward(u.val / v.val, (u.grad * v.val - v.grad * u.val) / v.val ^ 2);

end
end