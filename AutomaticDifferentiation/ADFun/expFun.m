function h = expFun(u)

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(exp(u), 0);
else
    h = ValGradForward(exp(u.val), exp(u.val)*u.Grad);
end

end