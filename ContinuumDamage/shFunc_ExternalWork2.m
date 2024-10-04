classdef shFunc_ExternalWork2 < handle
    
    properties (Access = public)

    end
    
    properties (Access = private)
        velocity 
        density 
        stress
        displacement 

    end
       
    methods (Access = public)
        
        function obj = shFunc_ExternalWork2(cParams)
            obj.
            
        end
        
        function energy = computeFunction(obj,quadOrder)
            
       
       
        end
        
        function jacobian = computeJacobian(obj,quadOrder)
            
          
            
        end
        
        function hessian = computeHessian(obj,quadOrder)
            
         
        
        end     
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
          
        end
        
    end
    
end