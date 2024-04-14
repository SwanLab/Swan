function solution=calc_ellipse(AOA,del_x,del_y)
options = optimoptions(@fsolve,'Algorithm','trust-region','OptimalityTolerance',1.0000e-8)
 
initial_guess=[0.001;0.001];
 
[solution] = fsolve(@(del_ab)ellipse_system(AOA,del_x,del_y,del_ab),initial_guess,options);


end