function h = minusFun(u,v)

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u - v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar

    [numRows,numCols] = size(u);

    if numRows > 1 && numCols > 1
        h = ValGradForward(u(1,:) - v.val(1,:),  - v.grad);

    else
        h = ValGradForward(u - v.val, - v.grad);

    end


elseif ~isa(v,'ValGradForward') %v is a scalar

    [numRows,numCols] = size(v);

    if numRows > 1 && numCols > 1
        h = ValGradForward(u.val(1,:) - v(1,:), u.grad);

    else
        h = ValGradForward(u.val - v, u.grad);

    end

else
    h = ValGradForward(u.val - v.val, u.grad - v.grad);

end
end