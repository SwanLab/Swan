classdef JacobiPreconditioner < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        D
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = JacobiPreconditioner(cParams)
            obj.init(cParams)
            
        end
        
        function x = apply(obj,r)
            x=obj.D\r;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.D = sparse(diag(diag(cParams.lhs)));
        end
        
    end
    
end