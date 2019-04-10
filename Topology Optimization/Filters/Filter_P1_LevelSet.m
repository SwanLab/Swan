classdef Filter_P1_LevelSet < Filter_P1 & Filter_LevelSet
    
    methods (Access = public)
        
        function obj = Filter_P1_LevelSet(cParams)
            obj@ Filter_P1(cParams);
            obj@ Filter_LevelSet(cParams);
        end
        
        function preProcess(obj)
            preProcess@Filter_P1(obj)
            preProcess@Filter_LevelSet(obj)
        end
        
    end
    
    methods (Access = protected)
        
        function x0 = computeP0fromP1(obj,x)
            M2 = obj.computeM2(x);
            x0 = obj.P_operator*M2;
        end
        
    end
    
    methods (Access = private)
        
        function M2 = computeM2(obj,x)
            switch obj.diffReacProb.geometry.type
                case 'TRIANGLE'
                    M2 = obj.faireF2(obj.mesh.coord',obj.mesh.connec',x);
                otherwise
                    M2 = obj.computeRHS(x,ones(size(x)));
            end
        end
        
    end
    
end
