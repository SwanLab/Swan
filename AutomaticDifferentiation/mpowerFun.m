function h = mpowerFun(u,v) %VALDER/MPOWER overloads ^ with at least one valder

if ~isa(u,'ValDerForward') %u is a scalar
    h = ValDerForward(u^v.val, u^v.val*log(u)*v.der);

elseif ~isa(v,'ValDerForward') %v is a scalar
    h = ValDerForward(u.val^v, v*u.val^(v-1)*u.der);

else
    h = exp(v*log(u)); %call overloaded log, * and exp

end
end