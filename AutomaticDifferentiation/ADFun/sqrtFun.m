function h = sqrtFun(u)

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(sqrt(u), 0);
else
    h = ValGradForward(sqrt(u.val), u.grad/(2*sqrt(u.val)));
end

end