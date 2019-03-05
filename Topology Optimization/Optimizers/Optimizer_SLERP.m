classdef Optimizer_SLERP < Optimizer_Unconstrained
    
    properties (Access = public)
        theta = 0.1
    end
    
    properties  (GetAccess = public, SetAccess = private)
        optimality_tol
    end
    
    properties (Access = private)
        normalizedPhi
        normalizedGrad
        coefPhi
        coefGrad
    end
    
    methods (Access = public)
        
        function obj = Optimizer_SLERP(settings,epsilon)
            obj@Optimizer_Unconstrained(settings,epsilon);
            obj.max_constr_change = +Inf;
            obj.nconstr = settings.nconstr;
        end
        
        function phi = computeX(obj,phi,g)
            obj.computeNormalizedLevelSet(phi);
            obj.computeNormalizedGradient(g);
            obj.computeTheta();
            obj.computeCoeficients();
            phi = obj.updateLevelSet();
            obj.updateOptimalityConditionValue();
        end
        
    end
    
    methods (Access = private)
        
        function computeNormalizedLevelSet(obj,phi)
            obj.normalizedPhi = obj.normalizeFunction(phi);
        end
        
        function computeNormalizedGradient(obj,g)
            obj.normalizedGrad = obj.normalizeFunction(g);
        end
        
        function computeCoeficients(obj)
            k = obj.line_search.kappa;
            t = obj.theta;
            obj.coefPhi  = sin((1-k)*t)/sin(t);
            obj.coefGrad = sin(k*t)/sin(t);
        end
        
        function phi = updateLevelSet(obj)
            cPhi = obj.coefPhi;
            cGra = obj.coefGrad;
            phiN = obj.normalizedPhi;
            g    = obj.normalizedGrad;
            phi = cPhi*phiN + cGra*g;
        end
        
        function computeTheta(obj)
            phiN = obj.normalizedPhi;
            g    = obj.normalizedGrad;
            phiXg = obj.scalar_product.computeSP(phiN,g);
            obj.theta = real(acos(phiXg));
        end
        
        function updateOptimalityConditionValue(obj)
            obj.opt_cond = obj.theta;
        end
        
        function x = normalizeFunction(obj,x)
            norm2 = obj.scalar_product.computeSP(x,x);
            xNorm = sqrt(norm2);
            x = x/xNorm;
        end
        
    end
    
    methods
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = (0.0175/1e-3)*obj.target_parameters.optimality_tol;
        end
        
    end
    
end