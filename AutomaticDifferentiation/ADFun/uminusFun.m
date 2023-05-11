function h = uminusFun(u)

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(u - v.val, v.Grad);

else
    h = ValGradForward(u.val - v.val, u.Grad - v.Grad);

end
end