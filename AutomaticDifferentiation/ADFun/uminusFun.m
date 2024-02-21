function h = uminusFun(u)

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(- u, 0);

else
    h = ValGradForward(- u.val, - u.grad);

end
end