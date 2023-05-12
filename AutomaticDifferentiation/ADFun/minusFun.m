function h = minusFun(u,v)

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u - v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(u - v.val, - v.grad);

elseif ~isa(v,'ValGradForward') %v is a scalar
    h = ValGradForward(u.val - v, u.grad);

else
    h = ValGradForward(u.val - v.val, u.grad - v.grad);

end
end