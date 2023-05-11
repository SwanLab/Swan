function h = mtimesFun(u,v) %ValGrad/MTIMES overloads * with at least one valGrad

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u * v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(u * v.val, u * v.grad);

elseif ~isa(v,'ValGradForward') %v is a scalar
    h = ValGradForward(v * u.val, v * u.grad);

else
    h = ValGradForward(u.val * v.val, u.grad * v.val + u.val * v.grad);

end
end