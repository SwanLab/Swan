classdef Optimizer_SLERP < Optimizer_Unconstrained
    
    properties (Access = public)
        theta = 0.1
    end
    
    properties (GetAccess = public, SetAccess = protected)
        type = 'SLERP'
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
        
        function obj = Optimizer_SLERP(settings)
            obj@Optimizer_Unconstrained(settings);
        end
        
        function compute(obj)
            phi = obj.designVariable.value;
            g   = obj.objectiveFunction.gradient;
            obj.computeNormalizedLevelSet(phi);
            obj.computeNormalizedGradient(g);
            obj.computeTheta();
            obj.computeCoeficients();
            phi = obj.updateLevelSet();
            obj.updateOptimalityConditionValue();
            obj.designVariable.value = phi;
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
            obj.theta = max(real(acos(phiXg)),1e-14);
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
    
    methods (Access = protected)
        
        function opt = obtainOptimalityTolerance(obj)
            opt = (0.0175/1e-3)*obj.targetParameters.optimality_tol;            
        end
        
    end
    
end