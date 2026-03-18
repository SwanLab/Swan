classdef CohesiveFunctional < handle
    
    properties (Access = public)
        mesh
        material
 
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CohesiveFunctional(cParams)
            obj.init(cParams)
            
        end



        function E = computeCost(obj,u,quadOrder)
           
        end

        function F = computeResidual(obj,u,quadOrder)
            
            obj.tractionLaw.
        end

        function H = computeDerivativeResidual(obj,u,quadOrder)
            %computeElementalCohesiveStiffnessMatrix()
            %assambeeMatrix()
            % ho fem amb la secant



            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.material = cParams.material;
        end

        function Felem = computeElementalCohesiveForceVector()
        end
        
    end
    
end