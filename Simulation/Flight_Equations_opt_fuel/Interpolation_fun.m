function [ALPHA_Ajusted, PI_Ajusted] = Interpolation_fun(ALPHA,PI,T_f,t,num_cons)
%Constants:
Sol_size = num_cons(1,3);

%Calculations:
T_solutions = linspace(0,T_f,Sol_size);

ALPHA_Ajusted = interp1(T_solutions,ALPHA,t);
PI_Ajusted = interp1(T_solutions,PI,t);

end

