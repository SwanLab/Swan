function h = mrdivideFun(u,v)

if ~isa(u,'ValGradForward') && ~isa(v,'ValGradForward') %u and v are scalars
    h = ValGradForward(u / v, 0);

elseif ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(u / v.val, (- v.grad * u) / v.val ^ 2);

elseif ~isa(v,'ValGradForward') %v is a scalar
    h = ValGradForward(u.val / v, u.grad / v);

else
    if ~isa(u.val,'ValGradForward') && ~isa(v.val,'ValGradForward') %u and v are scalars

        h = ValGradForward(u.val / v.val, (u.grad * v.val - v.grad .* u.val) / v.val ^ 2);

        % elseif ~isa(u,'ValGradForward') %u is a scalar
        %
        %     [numRows,numCols] = size(u);
        %
        %     if numRows > 1 && numCols > 1
        %         h = ValGradForward(v.val * u(1,1:end), v.grad .* u);

    else

        h = ValGradForward(u.val / v.val, (u.grad * v.val - v.grad * u.val) / v.val ^ 2);

    end

end
end