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



        function E = computeCost(obj,u,quadOrder, ????)
            

        end

        function F = computeGradient(obj,u,quadOrder, ????)
                
        end

        function H = computeHessian(obj,u,quadOrder,?)
        
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.material = cParams.material;
        end
        
    end
    
end