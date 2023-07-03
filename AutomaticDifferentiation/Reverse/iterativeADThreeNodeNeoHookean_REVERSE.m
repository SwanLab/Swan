function [val, grad] = iterativeADThreeNodeNeoHookean_REVERSE(u0)
    C1 = 1; %material constant
    D1 = 1; %material constant
    D = [1,-1,0; -1,2,-1; 0,-1,1];  
    [~,numElem] = size(u0); 
    G = zeros(1,numElem);
    G(1,numElem) = 1 * 10^(0);   
    u1 = ReverseVar(u0(1));
    u2 = ReverseVar(u0(2));
    u3 = ReverseVar(u0(3));
    
    F{1,1} = u1*D(1,1); 
    F{1,2} = u2*D(1,2); 
    F{1,3} = u3*D(1,3);
    F{2,1} = u1*D(1,2); 
    F{2,2} = u2*D(2,2); 
    F{2,3} = u3*D(2,3);
    F{3,1} = u1*D(1,3); 
    F{3,2} = u2*D(3,2); 
    F{3,3} = u3*D(3,3);
    
    I1 = F{1,1}^2 + F{2,2}^2+ F{3,3}^2;
    
    %J = F{1,1} * ( F{2,2} * F{3,3} - F{2,3} * F{3,2} ) - F{1,2} * ( F{2,1} * F{3,3} - F{2,3} * F{3,1} ) + F{1,3} * ( F{2,1} * F{3,2} - F{2,2} * F{3,1} );
    
    J = F{1,1} * F{2,2} * F{3,3};
    
    fun1 = ( I1 - 3 )*C1;
    
    fun2 =  log(J)*C1*2;
    
    fun3 = ( J - 1 )^2*D1;
    
    fun4 = u1*G(1)  + u2*G(2) + u3*G(3);
    
    fun = fun1 - fun2 + fun3 - fun4;
    
    val = fun.value;
    fun.grad_value = 1;
    grad = [u1.computeGradient(),u2.computeGradient(),u3.computeGradient()];
end