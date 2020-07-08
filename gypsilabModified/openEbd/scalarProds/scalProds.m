function[v,normH10_ab,scal01] = scalProds(a,kernel,rho)
% This function computes the vector which entries are
% v(i) = 2 \pi\int_{a(1)}^{a(2)}r(f'(r)ep'(r))
% where f' is the der property of kernel
% and e_p(r) = Cp(rho) J_0(rho(p) r)r
% Where Cp is a normalization constant in $H^1_0$. (See function 'Cp').
% The function uses oscillatory integrals through the matlab function
% "integral" (which is vectorized, but quite slow !).
% For big applications, we recommend using a custom Kernel object where the
% scalFunc is optimized (explicit if possible). 

warning off MATLAB:integral:NonFiniteValue


if nargin==0
    run('scalProdsVal.m')
    % Unitary Test
else
    scalFunc = kernel.scalFunc;
    normFunc = kernel.normFunc;
    
    if length(a) == 1
        b = 1;
    else
        b = a(2);
        a = a(1);
    end
    % Vector valued function
    % (i.e. when x is scalar, fun(x) is a vector (here of size
    % [1,length(rho)])
    v = scalFunc(a,b,rho);
    
    if a==0
        leftEdge = 0;
    else
        if kernel.singular
            leftEdge = Inf;
        else
            leftEdge = scalFunc(0,a,rho);
        end
    end
    
    if b==1
        rightEdge = 0;
    else
        rightEdge = scalFunc(b,1,rho);
    end
    
    scal01 = (leftEdge + v + rightEdge);
    normH10_ab = normFunc(a,b);
    
end

v = v(:);
scal01 = scal01(:);

warning on MATLAB:integral:NonFiniteValue

end