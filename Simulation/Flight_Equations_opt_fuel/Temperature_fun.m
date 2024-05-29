function Temp = Temperature_fun(x2,Earth_cons)
%Constants:
Temp_0 = Earth_cons(1,4);
Temp_11 = Earth_cons(1,5);
L_Temp = Earth_cons(1,6);

%Variable
h = x2;

%Calculations
if h < 11000
Temp = Temp_0-L_Temp.*h;

else
Temp = Temp_11;

end