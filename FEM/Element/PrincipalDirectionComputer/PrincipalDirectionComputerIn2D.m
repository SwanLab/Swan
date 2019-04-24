classdef PrincipalDirectionComputerIn2D < PrincipalDirectionComputer
    
    methods (Access = public)
        
        function obj = PrincipalDirectionComputerIn2D(cParams)
            obj.ndim = 2;
            obj.init();
        end
        
    end
    
    methods (Access = public)
        
        function compute(obj,tensor)
            s1  = squeeze(tensor(1,1,:));
            s2  = squeeze(tensor(1,2,:));
            s12 = squeeze(tensor(1,3,:));
            for i = 1:2
                for j = 1:2
                    obj.direction(i,j,:) = obj.directionFunction{i,j}(s1,s2,s12);
                end
            end
        end
        
    end
    
    methods (Access = protected)
           
        function obtainEigenVectors(obj)
            s1 = sym('s1','real');
            s2 = sym('s2','real');
            s12 = sym('s12','real');
            S = [s1 s12; s12 s2];
            [vS,~] = eig(S);
            obj.eigenVectors = vS;
        end        

    end
    
end