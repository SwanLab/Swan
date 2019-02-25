classdef Filter_P1_LevelSet < Filter_P1 & Filter_LevelSet
    
    methods (Access = public)
        
        function preProcess(obj)
            preProcess@Filter_P1(obj)
            preProcess@Filter_LevelSet(obj)
        end
        
        function x_gp = getP0fromP1(obj,x)
            if obj.xHasChanged(x)
                switch obj.diffReacProb.geometry.type
                    case 'TRIANGLE'
                        M2 = obj.faireF2(obj.mesh.coord',obj.mesh.connec',x);
                    otherwise
                        M2 = obj.computeRHS(x,ones(size(x)));
                end
                x_gp = obj.P_operator*M2;
                obj.x_reg = x_gp;
            else
                x_gp = obj.x_reg;
            end
        end
        
    end
    
    methods (Access = private)
        
        function itHas = xHasChanged(obj,x)
            itHas = ~norm(x) == norm(obj.x);
        end
    end
    
end
