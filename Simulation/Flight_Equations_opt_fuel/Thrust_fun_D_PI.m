function Thrust_D_PI = Thrust_fun_D_PI(rho_eng,PI_Ajusted,Engine_cons,Earth_cons)
    %Constants:
    T_max = Engine_cons(1,1);
    T_idle = Engine_cons(1,2);
    rho_0 = Earth_cons(1,2);

    Thrust_D_PI = (rho_eng./rho_0).*(T_max-T_idle);
end

