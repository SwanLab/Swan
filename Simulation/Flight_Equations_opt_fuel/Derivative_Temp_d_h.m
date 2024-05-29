function D_Temp_D_h = Derivative_Temp_d_h(x2,Earth_cons)
 %Constants:
 L_Temp = Earth_cons(1,6);
 h = x2;

 %Calculations
 if h<11000
 D_Temp_D_h = -L_Temp;
 else
 D_Temp_D_h = 0;
 end
end

