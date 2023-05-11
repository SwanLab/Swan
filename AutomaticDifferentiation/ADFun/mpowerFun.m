function h = mpowerFun(u,v) %VALGrad/MPOWER overloads ^ with at least one valGrad

if ~isa(u,'ValGradForward') %u is a scalar
    h = ValGradForward(u^v.val, u^v.val*log(u)*v.grad);

elseif ~isa(v,'ValGradForward') %v is a scalar
    h = ValGradForward(u.val^v, v*u.val^(v-1)*u.grad);

else
    h = exp(v*log(u)); %call overloaded log, * and exp

end
end