function h = mtimesFun(u,v) %ValGrad/MTIMES overloads * with at least one valGrad

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u * v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar

    [numRows,numCols] = size(u);

    if numRows > 1 && numCols > 1
        h = ValGradForward(v.val * u(1,1:end), v.grad .* u);

    else
        h = ValGradForward(u * v.val, u * v.grad);

    end

elseif ~isa(v,'ValGradForward') %v is a scalar

    [numRows,numCols] = size(v);

    if numRows > 1 && numCols > 1
        h = ValGradForward(u.val * v(1,1:end), u.grad .* v);

    else
        h = ValGradForward(v * u.val, v * u.grad);

    end

else
    if ~isa(u.val,'ValGradForward') && ~isa(v.val,'ValGradForward') %u and v are scalars
        
        h = ValGradForward(u.val * v.val, u.grad .* v.val + u.val * v.grad);

        % elseif ~isa(u,'ValGradForward') %u is a scalar
        %
        %     [numRows,numCols] = size(u);
        %
        %     if numRows > 1 && numCols > 1
        %         h = ValGradForward(v.val * u(1,1:end), v.grad .* u);

    else
        
        h = ValGradForward(u.val * v.val, u.grad * v.val + u.val * v.grad);

    end

end
end