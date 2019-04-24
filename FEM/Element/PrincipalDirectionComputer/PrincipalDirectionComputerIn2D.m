classdef PrincipalDirectionComputerIn2D < PrincipalDirectionComputer
    
    methods (Access = public, Static)
        
        function obj = PrincipalDirectionComputerIn2D(cParams)
            obj.init();
        end
        
    end
    
    methods (Access = public)
        
        function init(obj)
            s1 = sym('s1','real');
            s2 = sym('s2','real');
            s12 = sym('s12','real');            
            S = [s1 s12; s12 s2];
            [~,vS] = eig(S);
            for i = 1:2
                for j = 1:2
                 obj.directionFunction{i,j} = obj.symbolic2Function(vS(i,j));
                end
            end            
        end
        
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
    
    methods (Access = private, Static)
        
        function sF = symbolic2Function(s)
            if isequal(s,sym(0))
                sF = @(s1,s2,s12) 0;
            else
                sF = matlabFunction(s);
            end
        end
    end
    
end