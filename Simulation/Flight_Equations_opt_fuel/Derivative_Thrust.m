function [Thrust_D_h, Thrust_D_v] = Derivative_Thrust(PI_Ajusted,Engine_cons,Earth_cons,D_rho_eng_D_v,D_rho_eng_D_h)
    %Constants:
    T_max = Engine_cons(1,1);
    T_idle = Engine_cons(1,2);
    rho_0 = Earth_cons(1,2);

    %Calculations:

    Thrust_0 = (T_max-T_idle).*PI_Ajusted + T_idle;

    %Derivative of speed:
    a = Thrust_0./rho_0;

    Thrust_D_v = a.*D_rho_eng_D_v;

    %Derivative of height:
    Thrust_D_h = a.*D_rho_eng_D_h;
end

