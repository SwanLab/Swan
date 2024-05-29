function [D_rho_eng_D_h, D_rho_eng_D_v] = Derivative_airdensity_eng_d_h(M,rho,Earth_cons,D_M_D_v,D_M_D_h,D_rho_D_h)
    %Constants:
    gamma = Earth_cons(1,8); %Gamma of air
    
    %Calculations:
    %Derivation of rho_eng relative to speed:
    Base = 1 +0.5.*(gamma-1).*M.^2;
    a = (1./(gamma-1));
    b = (gamma-1).*M;
    Exponent = (1./(gamma-1))-1;
    D_rho_eng_D_v = rho.*a.*b.*Base.^(Exponent).*D_M_D_v;

    
    %Derivation of rho_eng relative to height:
    First_term = Base.^(a);
    Second_term = rho.*a.*b.*Base.^(Exponent);

    D_rho_eng_D_h = First_term.*D_rho_D_h + Second_term.*D_M_D_h;
end

