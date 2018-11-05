classdef Filter_P1_LevelSet < Filter_P1
    methods (Abstract)
        preProcess(obj)
    end
    
    methods
        function x_gp = getP0fromP1(obj,x)
            xHasChanged = ~norm(x) == norm(obj.x);
            
            if xHasChanged
                switch obj.diffReacProb.geometry.type
                    case 'TRIANGLE'
                        M2 = obj.faireF2(obj.mesh.coord',obj.mesh.connec',x);
                    otherwise
                        obj.unfitted_mesh.computeMesh(x);
                        M2 = obj.computeRHS(ones(size(x)));
                end
                x_gp = obj.P_operator*M2;
                obj.x_reg = x_gp;
            else
                x_gp = obj.x_reg;
            end
        end
    end
end