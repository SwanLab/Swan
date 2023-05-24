function h = minusFun(u,v)

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u - v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar
    
    if isequal(size(u),[3,3])
        h = ValGradForward(u(1,:) - v.val(1,:),  - v.grad);

    else
        h = ValGradForward(u - v.val, - v.grad);

    end


elseif ~isa(v,'ValGradForward') %v is a scalar

    if isequal(size(v),[3,3])
        h = ValGradForward(u.val(1,:) - v(1,:), u.grad);

    else
        h = ValGradForward(u.val - v, u.grad);

    end

else
        h = ValGradForward(u.val - v.val, u.grad - v.grad);

end
end