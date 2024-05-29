function [D_CL_D_h, D_CL_D_v] = Derivative_CL(D_CL_alpha_D_h,D_CL_alpha_D_v,ALPHA_Ajusted,CL_alpha)
    %Calculations:
    %Relative to speed:
    D_CL_D_v = ALPHA_Ajusted.*D_CL_alpha_D_v;

    %Relative to height:
    D_CL_D_h = ALPHA_Ajusted.*D_CL_alpha_D_h;

end

