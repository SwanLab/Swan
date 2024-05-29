function Thrust = Thrust_fun(rho_eng,PI_Ajusted,Engine_cons,Earth_cons)
%Constants:
T_max = Engine_cons(1,1);
T_idle = Engine_cons(1,2);
rho_0 = Earth_cons(1,2);

%Calculations:

Thrust_0 = (T_max-T_idle).*PI_Ajusted + T_idle;

Thrust = (rho_eng./rho_0).*Thrust_0;

end

