classdef MinimumEigenValueFunctional < handle
    
    properties (Access = public)
        value
        gradient
    end
    
    properties (Access = private)
        density
    end
    
    properties (Access = private)
       eigModes 
       designVariable
    end
    
    methods (Access = public)
        
        function obj = MinimumEigenValueFunctional(cParams)
            obj.init(cParams)
        end
        
        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            % obj.designVariable.update(x);
            obj.computeDensity(); 
            [f,dfdx]= obj.eigModes.computeFunctionAndGradient(obj.density);
            obj.value = f;  
            obj.gradient = dfdx;            
        end
        
        function computeFunction(obj)
            obj.computeDensity();
            [f,~]= obj.eigModes.computeFunctionAndGradient(obj.density);
            obj.value = f;
            %obj.normalizeFunction();
        end
        
        function computeGradient(obj)
            obj.computeDensity();
            [~,dfdx] = obj.eigModes.computeFunctionAndGradient(obj.density);
            obj.gradient = dfdx;
            %obj.normalizeGradient();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
        end
        
        function computeDensity(obj)
            densDomain  = obj.designVariable.fun;
            s.operation = @(xV) obj.computeComplementaryDensity(densDomain,xV);
            densHole = DomainFunction(s);
            obj.density = densHole;
        end
        
        function rhoP = computeComplementaryDensity(obj,fun,xV)
            rho = 1 - fun;
            s.beta = 2.0;
            s.eta  = 0.5;
            projector = HeavisideProjector(s);
            rhoPVal = projector.project(rho);
            rhoP = LagrangianFunction.create(obj.designVariable.fun.mesh,1,fun.order);
            rhoP.fValues = rhoPVal;

            rho3 = round(rho.fValues);
            rho3 = max(0,min(1,rho3));
            rhoP3 = LagrangianFunction.create(obj.designVariable.fun.mesh,1,fun.order);
            rhoP3.fValues = rho3;
            
            s.beta = 2.0;
            s.eta = 0.05;
            s.filterStep = 'LUMP';
            s.mesh  = obj.designVariable.fun.mesh;
            s.trial = LagrangianFunction.create(obj.designVariable.fun.mesh,1,'P1');
            filter = FilterAndProject(s);
            rhoP4 = filter.compute(rho, 3);
         end

%         function createHeavisideProjector(obj, beta, eta)
%             obj.HeavisideProjector()

    end
    
    
end