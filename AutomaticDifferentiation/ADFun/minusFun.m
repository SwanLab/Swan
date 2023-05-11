function h = minusFun(u,v)

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u - v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(u - v.val, v.Grad);

elseif ~isa(v,'ValGradForward') %v is a scalar
    h = ValGradForward(v - u.val, u.Grad);

else
    h = ValGradForward(u.val - v.val, u.Grad - v.Grad);

end
end