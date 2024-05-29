function [D_M_D_h, D_M_D_v] = Derivative_Mach(v,D_Temp_D_h,Temp,Earth_cons)
   %Constants:
    r_gas= Earth_cons(1,7);
    gamma_air = Earth_cons(1,8);

    %Calculations:
    %Derivation of M relative to speed:
    Denominator = sqrt(r_gas.*gamma_air.*Temp);
    D_M_D_v = 1./Denominator;
    
    %Derivation of M relative to height:
    a =- 0.5.*v.*gamma_air.*r_gas;
    b = gamma_air.*r_gas.*Temp;
    D_M_D_h = a.*D_Temp_D_h.*b^(-1.5);
end

