function h = mtimesFun(u,v) %ValGrad/MTIMES overloads * with at least one valGrad

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u * v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar

    if isequal(size(u),[3,3])
        h = ValGradForward(v.val * u(1,1:3), v.grad .* u);

    else
        h = ValGradForward(u * v.val, u * v.grad);

    end

elseif ~isa(v,'ValGradForward') %v is a scalar

    if isequal(size(v),[3,3])
        h = ValGradForward(u.val * v(1,1:3), u.grad .* v);

    else
        h = ValGradForward(v * u.val, v * u.grad);

    end

else
    h = ValGradForward(u.val * v.val, u.grad * v.val + u.val * v.grad);

end
end