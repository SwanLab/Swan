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
            L = obj.Lchol;
            z = L\r;
            z = (L')\z;
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.LHS = cParams.LHS;            
        end
        
        function computeCholeskyFactorization(obj)
            opts.type     = 'ict';
            opts.droptol  = 1e-3;
            opts.diagcomp = 1e-2;
            L = ichol(obj.LHS, opts);        
            obj.Lchol = L ;
        end
    end
    
end