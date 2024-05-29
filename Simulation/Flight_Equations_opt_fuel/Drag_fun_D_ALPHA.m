function Drag_D_ALPHA = Drag_fun_D_ALPHA(v,rho,ALPHA_Ajusted,CL_alpha,Aircraft_cons,Aero_cons)
    %Constants:
    S = Aircraft_cons(1,1);
    K = Aero_cons(1,1);

    Drag_D_ALPHA = ALPHA_Ajusted.*K.*S.*rho.*CL_alpha.*v.^2;
end

