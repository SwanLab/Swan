classdef PrincipalStressDirections < handle
    
    properties (GetAccess = public, SetAccess = private)
        principalStressDir
        principalStress
    end
    
    properties (Access = private)
        stress
    end
    
    methods (Access = public)
        
        function obj = PrincipalStressDirections(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeStressBase();
            obj.normalizePrincipalDirection();
            obj.transformBaseToElementalBase();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.stress = cParams.stress;
        end
        
        function computeStressBase(obj)
            s = obj.stress;
            S = [s(1) s(3);s(3) s(2)];
            [V,D] = eig(S);
            obj.principalStressDir = V;
            obj.principalStress = D;
        end
        
    end
    
    methods (Access = private)
        
        function normalizePrincipalDirection(obj)
            V = obj.principalStressDir;
            V(:,1) = V(:,1)/norm(V(:,1));
            V(:,2) = V(:,2)/norm(V(:,2));
            obj.principalStressDir = V;
        end
        
        function transformBaseToElementalBase(obj)
            V = obj.principalStressDir;
            V(:,2) = V(:,2)*det(V);
            obj.principalStressDir = V;
        end
        
    end
end

