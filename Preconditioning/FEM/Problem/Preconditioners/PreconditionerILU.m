classdef PreconditionerILU < handle
    
    properties (Access = private)
        Lchol
    end
    
    properties (Access = private)        
        LHS
    end
    
    methods (Access = public)
        
        function obj = PreconditionerILU(cParams)
            obj.init(cParams)
            obj.computeCholeskyFactorization()
        end
        
        function z = apply(obj,r)
%             L = obj.Lchol;
%             z = L\r;
%             z = (L')\z;
            z = obj.Lchol' \ (obj.Lchol \ r);
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.LHS = cParams.LHS;            
        end
        
        function computeCholeskyFactorization(obj)
            L = ichol(obj.LHS);        
            obj.Lchol = L ;
        end
    end
    
end