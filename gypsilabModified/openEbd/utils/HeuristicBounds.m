function[doup] = HeuristicBounds(a,tol)

if nargin==0
    HeuristicBoundsVal();
else
    % Empirical values !
    gamma1_up = 0.31; gamma1_do = 0.28;
    gamma2_up = -0.14; gamma2_do = -0.6;
    
    gamma_up = max(gamma2_up+gamma1_up*log(1/tol),0.5);
    up = fix(gamma_up*1/a)+1;
    do = max(fix((gamma2_do+gamma1_do*log(1/tol))*1/a)-1,1);
    doup = [do,up];
end




end