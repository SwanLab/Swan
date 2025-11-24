classdef PreconditionerGaussSeidel < handle
    
    properties (Access = private)
        L 
        U 
        D
    end
    
    properties (Access = private)        
        LHS
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = PreconditionerGaussSeidel(cParams)
            obj.init(cParams)
            obj.computeLUFactorization()
        end
        
        function z = apply(obj,r)            
            Pl   = @(r) (obj.L + obj.D)\r;
            Pu   = @(r) (obj.U + obj.D)\r;
            LHSf = @(r) obj.LHS*r;
            z    = Preconditioner.multiplePrec(r,Pl,Pu,LHSf);            
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.LHS = cParams.LHS;            
        end
        
        function computeLUFactorization(obj)
            obj.L = tril(obj.LHS);
            obj.U = triu(obj.LHS);
            obj.D = diag(diag(obj.LHS));
        end
        
        
    end
    
    
end