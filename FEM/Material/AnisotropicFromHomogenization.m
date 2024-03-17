classdef AnisotropicFromHomogenization < Material
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        vadVariables
        Ctensor 
        sMesh
    end
    
    properties (Access = private)
        microParams
        fileName
    end
    
    methods (Access = public)
        
        function obj = AnisotropicFromHomogenization(cParams)
            obj.init(cParams)
   
        end
        
        function C = evaluate(obj,xV)
            C = obj.computeValues(xV);            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.microParams = cParams.microParams;
           obj.fileName    = cParams.fileName;           
        end
        
 

                 
        
    end
    
end