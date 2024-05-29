function Lift_D_ALPHA = Lift_fun_D_ALPHA(v,rho,ALPHA_Ajusted,CL_alpha,Aircraft_cons)
    %Constants:
    S = Aircraft_cons(1,1);

    Lift_D_ALPHA = 0.5.*S.*rho.*CL_alpha.*v.^2;
end

