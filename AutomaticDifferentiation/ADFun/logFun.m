function h = logFun(u)

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(log(u), 0);
else
    h = ValGradForward(log(u.val), u.der/u.val);
end

end