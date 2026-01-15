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
            % opts.type = 'ict';
            % opts.droptol = 1e-3;
            % opts.diagcomp = 1e-3;   % or larger if needed
            % L = ichol(obj.LHS, opts);
            % setup.type = 'ilutp';
            % setup.droptol = 1e-3;
            % [L,U] = ilu(obj.LHS,setup);

            opts.type     = 'ict';
            opts.droptol  = 1e-3;
            opts.diagcomp = 1e-2;   % try 1e-6 → 1e-2 if needed
            L = ichol(obj.LHS, opts);
            % L = ichol(obj.LHS);
            obj.Lchol = L ;
        end
    end
    
end