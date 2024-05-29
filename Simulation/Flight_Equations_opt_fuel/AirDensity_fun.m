function rho = AirDensity_fun(x2,Earth_cons)
%Constants:
g = Earth_cons(1,1);
rho_0 = Earth_cons(1,2);
rho_11 = Earth_cons(1,3);
Temp_0 = Earth_cons(1,4);
Temp_11 = Earth_cons(1,5);
L_Temp = Earth_cons(1,6);
r_gas= Earth_cons(1,7);


%Variable
h = x2;

%Calculations
if h < 11000
Base = 1-(L_Temp./Temp_0).*h;
Exponent = (g./(r_gas.*L_Temp))-1;

rho = rho_0.*Base.^Exponent;
else
Exponent = -(g.*(h-11000))./(r_gas.*Temp_11);

rho = rho_11.*exp(Exponent);
end