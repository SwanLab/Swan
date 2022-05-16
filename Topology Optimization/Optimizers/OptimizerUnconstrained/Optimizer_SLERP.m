classdef Optimizer_SLERP < Optimizer_Unconstrained
    
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
        theta
    end
    
    methods (Access = public)
        
        function obj = Optimizer_SLERP(cParams)
            obj@Optimizer_Unconstrained(cParams);
        end
        
        function compute(obj)
            obj.computeNormalizedLevelSet();
            obj.computeNormalizedGradient();
            obj.computeTheta();
            obj.computeCoeficients();
            phi = obj.updateLevelSet();
            obj.updateOptimalityConditionValue();
            obj.designVariable.update(phi);
        end
        
        function storeConvergenceVariablesValues(obj)
            obj.convergenceVars.reset();
            obj.convergenceVars.append(obj.incF);
            obj.convergenceVars.append(obj.incX);
            obj.convergenceVars.append(obj.lineSearch.value);
            obj.convergenceVars.append(obj.lineSearch.nTrials);
            obj.convergenceVars.append(obj.obtainThetaValue());
        end
        
    end
    
    methods (Access = private)
        
        function computeNormalizedLevelSet(obj)
            phi = obj.designVariable.value;
            obj.normalizedPhi = obj.normalizeFunction(phi);
        end
        
        function computeNormalizedGradient(obj)
            g   = obj.objectiveFunction.gradient;
            obj.normalizedGrad = obj.normalizeFunction(g);
        end
        
        function computeCoeficients(obj)
            k = obj.lineSearch.value;
            t = obj.theta;
            obj.coefPhi  = sin((1-k)*t)/sin(t);
            obj.coefGrad = sin(k*t)/sin(t);
        end
        
        function phi = updateLevelSet(obj)
            cPhi = obj.coefPhi;
            cGra = obj.coefGrad;
            phiN = obj.normalizedPhi;
            g    = obj.normalizedGrad;
            phi  = cPhi*phiN + cGra*g;
        end
        
        function computeTheta(obj)
            phiN      = obj.normalizedPhi;
            g         = obj.normalizedGrad;
            phiXg     = obj.scalar_product.computeSP(phiN,g);
            obj.theta = max(real(acos(phiXg)),1e-14);
        end
        
        function t = obtainThetaValue(obj)
            if isempty(obj.theta)
                obj.computeNormalizedLevelSet();
                obj.computeNormalizedGradient();
                obj.computeTheta(); 
            end
            t = obj.theta*180/pi;
        end
        
        function updateOptimalityConditionValue(obj)
            obj.optimalityCond = obj.thetaToNorm(obj.theta);
        end
        
        function x = normalizeFunction(obj,x)
            norm2 = obj.scalar_product.computeSP(x,x);
            xNorm = sqrt(norm2);
            x = x/xNorm;
        end
        
    end
    
    methods (Access = private, Static)
        
        function norm = thetaToNorm(theta)
            y0 = 1e-4;
            y1 = 1e-3;
            x0 = 0.1*pi/180;
            x1 = 1*pi/180;
            f = @(x) y0 + (y1-y0)/(x1-x0)*(x-x0);
            norm = f(theta);
        end
        
    end
    
    
end