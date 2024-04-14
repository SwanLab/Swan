classdef calc_ellipse_classe
    properties
        AOA
        del_x
        del_y
    end
    
    methods 
        function obj = calc_ellipse_classe(AOA,del_x,del_y) %Funció constructora que assigna els paràmetres d'entrada a cada propietat de la classe
            obj.AOA = AOA;
            obj.del_x = del_x;
            obj.del_y = del_y;
        end
        
        function solution=solvesys(obj)
           options = optimoptions(@fsolve,'Algorithm','trust-region','OptimalityTolerance',1.0000e-8);
 
           initial_guess=[0.001;0.001];
 
           [solution] = fsolve(@(del_ab)ellipse_system(obj.AOA,obj.del_x,obj.del_y,del_ab),initial_guess,options);
        end


        function F=ellipse_system(AOA,del_x,del_y,del_ab)

            F=[cos(atan(del_ab(2)/del_ab(1))+AOA)*sqrt((del_ab(1))^2+(del_ab(2))^2)-del_x;
               sin(atan(del_ab(2)/del_ab(1))+AOA)*sqrt((del_ab(1)^2)+(del_ab(2)^2))-del_y];
        
        end


    end
end

