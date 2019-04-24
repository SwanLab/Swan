classdef PrincipalDirectionComputerIn3D < PrincipalDirectionComputer
    
    methods (Access = public, Static)
        
        function obj = PrincipalDirectionComputerIn3D(cParams)
            obj.init()
        end
        
    end
    
    methods (Access = public)
        
        function init(obj)
            s1 = sym('s1','real');
            s2 = sym('s2','real');
            s3 = sym('s3','real');
            s12 = sym('s12','real');
            s13 = sym('s13','real');
            s23 = sym('s23','real');
            S = [s1 s12 s13; s12 s2 s23; s13 s23 s3];
            [~,vS] = eig(S);            
            for i = 1:3
                for j = 1:3
                 obj.directionFunction{i,j} = obj.symbolic2Function(vS(i,j));
                end
            end            
        end
        
        function compute(obj,tensor)
            s1  = squeeze(tensor(1,1,:));
            s2  = squeeze(tensor(1,2,:));
            s3  = squeeze(tensor(1,3,:));
            s12 = squeeze(tensor(1,4,:));
            s13 = squeeze(tensor(1,5,:));
            s23 = squeeze(tensor(1,6,:));            
            for i = 1:3
                for j = 1:3
                 d = obj.directionFunction{i,j}(s1,s2,s3,s12,s13,s23);
                 obj.direction(i,j,:) = round(d,14);
                end
            end
        end
        
    end    
    
    methods (Access = private, Static)
        
        function sF = symbolic2Function(s)
            if isequal(s,sym(0))
                sF = @(s1,s2,s3,s12,s13,s23) 0;
            else
                sF = matlabFunction(s);
            end
        end
    end    
        
end