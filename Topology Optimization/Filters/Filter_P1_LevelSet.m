classdef Filter_P1_LevelSet < Filter_P1
    methods (Abstract)
        preProcess(obj)
    end
    
    methods
        function obj = Filter_P1_LevelSet(problemID,scale)            
            obj@Filter_P1(problemID,scale);
        end        
        
        function x_gp = getP0fromP1(obj,x)
            if norm(x) == norm(obj.x)
                x_gp=obj.x_reg;
            else
                switch obj.geometry.type
                    case 'TRIANGLE'
                        M2=obj.faireF2(obj.coordinates',obj.connectivities',x);
                    otherwise
                        M2=obj.computeRHS(x);
                end
                x_gp = obj.P_operator*M2;
                obj.x_reg=x_gp;
            end
        end
    end
end

