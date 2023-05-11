function h = cosFun(u)

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(cos(u), 0);
else
    h = ValGradForward(cos(u.val), -sin(u.val)*u.Grad);
end

end