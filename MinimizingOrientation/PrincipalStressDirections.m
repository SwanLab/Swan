classdef PrincipalStressDirections < handle
    
    properties (GetAccess = public, SetAccess = private)
        principalStressDir   
        stress
    end
    
    properties (Access = private)
        strain
        Ctensor
    end
    
    methods (Access = public)
        
        function obj = PrincipalStressDirections(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeStress();
            obj.computeStressPrincipalDirections();                        
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.strain = cParams.strain;
            obj.Ctensor = cParams.Ctensor;            
        end
        
        
        function computeStress(obj)
            C = obj.Ctensor;
            e = obj.strain;
            obj.stress = C*e;
        end
        
        function computeStressPrincipalDirections(obj)
            s = obj.stress;
            S = [s(1) s(3);s(3) s(2)];
            [V,D] = eig(S);
            V(:,1) = V(:,1)/norm(V(:,1));
            V(:,2) = V(:,2)/norm(V(:,2));
            obj.principalStressDir = V;
        end
        
        
    end
end

