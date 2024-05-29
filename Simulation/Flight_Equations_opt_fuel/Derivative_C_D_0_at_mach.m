function [D_C_D_0_D_h, D_C_D_0_D_v] = Derivative_C_D_0_at_mach(M,Aero_cons,D_M_D_v,D_M_D_h)
    %Constants:
    b = Aero_cons(1,5);
    c = Aero_cons(1,6);
    d = Aero_cons(1,7);
    e = Aero_cons(1,8);

    %Previous calculations:

    D_C_D_0_M = b + 2.*c.*M +3.*d.*M.^2 + 4.*e.*M.^3;

    %Derivative of CD_0 vs speed:

    D_C_D_0_D_v = D_C_D_0_M.*D_M_D_v;

    %Derivative of CD_0 vs height:

    D_C_D_0_D_h = D_C_D_0_M.*D_M_D_h;    

end

