function D_rho_D_h = Derivative_airdensity_d_h(x2,Earth_cons)
    %Constants:
    g = Earth_cons(1,1);
    rho_0 = Earth_cons(1,2);
    rho_11 = Earth_cons(1,3);
    Temp_0 = Earth_cons(1,4);
    Temp_11 = Earth_cons(1,5);
    L_Temp = Earth_cons(1,6);
    r_gas= Earth_cons(1,7);
    h = x2;

    if h<11000
        a = g./(r_gas.*L_Temp)-1;
        b = 1-(L_Temp./Temp_0).*h;
        c = (g./(r_gas.*L_Temp))-2;
        d = (-L_Temp./Temp_0);
    
        D_rho_D_h = rho_0.*a.*d.*b.^c;
    else
        Exponent = -(g.*(h-11000))./(r_gas.*Temp_11);
        b = -g./(r_gas.*Temp_11);
        D_rho_D_h = rho_11.*b.*exp(Exponent);
    end
end

