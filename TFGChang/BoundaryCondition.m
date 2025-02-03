classdef BoundaryCondition < handle
    
    properties (Access = public)
        forcesFormula
        dirConditions
        dirDofs
        nodesConditions
    end
    
    properties (Access = private)
        mesh
        uMesh 
        velocityFun
        pressureFun
    end

    methods (Access = public)
        
        function obj = BoundaryCondition(cParams)
            obj.init(cParams);
        end

        function compute(obj)

            obj.forcesFormula
            obj.dirConditions
            obj.dirDofs
            obj.nodesConditions
            %falta per completar

        end
        
    end
        
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.uMesh        = cParams.uMesh;
            obj.velocityFun  = cParams.velocityFun;
            obj.pressureFun  = cParams.pressureFun;
        end
    end

end