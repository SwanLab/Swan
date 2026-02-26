classdef PreconditionerJacobi < handle
    
    properties (Access = private)
        D
    end
    
    properties (Access = private)        
        LHS
    end
    
    methods (Access = public)
        
        function obj = PreconditionerJacobi(cParams)
            obj.init(cParams)
            obj.computeDiagonalMatrix()
        end
        
        function z = apply(obj,r)            
            Pl = @(r) (obj.D)\r;
            z  = Pl(r);
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.LHS = cParams.LHS;            
        end
        
        function computeDiagonalMatrix(obj)
            obj.D = diag(diag(obj.LHS));
        end
        
    end
    
end