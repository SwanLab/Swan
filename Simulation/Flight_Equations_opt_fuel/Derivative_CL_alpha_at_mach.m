function [D_CL_alpha_D_h, D_CL_alpha_D_v] = Derivative_CL_alpha_at_mach(M,Aero_cons,Aircraft_cons,D_M_D_v,D_M_D_h)
    %Constants:
    Cl_alpha_0 = Aero_cons(1,2);
    e_factor = Aero_cons(1,3);
    AR = Aircraft_cons(1,3);

    %Previous Calculations:

    Denominator = 1 + (Cl_alpha_0)./(e_factor.*AR.*pi);

    CL_alpha_0 = Cl_alpha_0./Denominator;

    %Derivation CL_Alpha vs Mach:

    Base = 1-M.^2;
    Exponent = -1.5;
    D_CL_alpha_D_M = CL_alpha_0.*M.*Base.^Exponent;

    %Derivation CL_Alpha vs speed:
    D_CL_alpha_D_v = D_CL_alpha_D_M.*D_M_D_v;

    %Derivation CL_Alpha vs height:
    D_CL_alpha_D_h = D_CL_alpha_D_M.*D_M_D_h;
end

