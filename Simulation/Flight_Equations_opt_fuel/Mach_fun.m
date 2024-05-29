function M = Mach_fun(Temp,v,Earth_cons)
%Constants:
r_gas= Earth_cons(1,7);
gamma_air = Earth_cons(1,8);

%Calculations
Denominator = sqrt(gamma_air.*r_gas.*Temp);

M = v./Denominator;
end