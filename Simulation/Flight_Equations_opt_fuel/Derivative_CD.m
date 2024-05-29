function [D_CD_D_h, D_CD_D_v] = Derivative_CD(D_CL_alpha_D_h,D_CL_alpha_D_v,D_C_D_0_D_h,D_C_D_0_D_v,ALPHA_Ajusted,CL_alpha,Aero_cons)
    %Constants:
    K = Aero_cons(1,1);
    
    %Derivation of speed:
    First_term = D_C_D_0_D_v;
    Second_term = 2.*K.*CL_alpha.*ALPHA_Ajusted.^(2).*D_CL_alpha_D_v;

    D_CD_D_v = First_term+Second_term;

    %Derivation of height:
    First_term = D_C_D_0_D_h;
    Second_term = 2.*K.*CL_alpha.*ALPHA_Ajusted.^(2).*D_CL_alpha_D_h;

    D_CD_D_h = First_term+Second_term;
end

