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

           %  dens.fValues = round(dens.fValues);  
           % dens.fValues = max(0,min(1,dens.fValues));
            obj.density = densHole;
        end
        
        function rho = computeComplementaryDensity(obj,fun,xV)
            rho = fun.evaluate(xV);
            rho = 1 - rho;
         end

    end
    
    
end