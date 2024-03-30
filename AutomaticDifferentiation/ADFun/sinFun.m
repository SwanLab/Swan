function h = sinFun(u)

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(sin(u), 0);
else
    h = ValGradForward(sin(u.val), cos(u.val)*u.grad);
end

end