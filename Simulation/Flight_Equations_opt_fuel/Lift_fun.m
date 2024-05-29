function Lift  = Lift_fun(C_L,rho,v,Aircraft_cons)
%Constants:
S = Aircraft_cons(1,1);

Lift = 0.5.*rho.*S.*C_L.*v.^2;

end

