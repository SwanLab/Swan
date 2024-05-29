function [Lift_D_h, Lift_D_v] = Derivative_Lift(v,rho,C_L,Aircraft_cons,D_CL_D_v,D_CL_D_h,D_rho_D_h)
    %Constants:
    S = Aircraft_cons(1,1);

    %Derivative with speed:
    First_term = C_L.*S.*rho.*v;
    Second_term =  0.5.*S.*rho.*v.^(2).*D_CL_D_v;

    Lift_D_v = First_term + Second_term;

    %Derivative with height:
    First_term = 0.5.*C_L.*S.*v.^(2).*D_rho_D_h;
    Second_term = 0.5.*S.*rho.*v.^(2).*D_CL_D_h;

    Lift_D_h = First_term + Second_term;
end

