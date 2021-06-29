classdef ScalarProduct < handle
    
    properties (Access = public)
        epsilon
        Ksmooth
        Msmooth
    end
    
    properties (Access = private)
       nVariables 
    end
    
    methods (Access = public)
        
        function obj = ScalarProduct(cParams)
            obj.init(cParams);
            obj.createMatrices(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function sp = computeSP(obj,f,g)
            spM = obj.computeSP_M(f,g);
            spK = obj.computeSP_K(f,g);
            sp  = obj.epsilon^2*spK + spM;
        end
        
        function sp = computeSP_M(obj,f,g)
            sp = obj.computeProduct(obj.Msmooth,f,g);
        end
        
        function sp = computeSP_K(obj,f,g)
            sp = obj.computeProduct(obj.Ksmooth,f,g);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.epsilon = cParams.epsilon;
            obj.nVariables = cParams.nVariables;
        end
        
        function createMatrices(obj,cParams)
            s = cParams.femSettings;
            s.mesh = cParams.mesh;
            physProb = DiffReact_Problem(s);
            physProb.preProcess();
            obj.Ksmooth = physProb.element.K;
            obj.Msmooth = physProb.element.M;
        end
        
        function n = computeProduct(obj,K,f,g)
            nx = length(f)/obj.nVariables;
            n = 0;
            for ivar = 1:obj.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                fs = f(i0:iF);
                gs = g(i0:iF);                
                n = n + fs'*K*gs;
            end              
        end
        
    end    
    
end
