classdef MinimumEigenValueFunctional < handle
    
    properties (Access = public)
        value
        gradient
    end
       
    properties (Access = private)
       density
       eigModes 
       designVariable
       mesh
       filter
       iter
    end
    
    methods (Access = public)
        
        function obj = MinimumEigenValueFunctional(cParams)
            obj.init(cParams)
            obj.createFilterAndProject()
        end
        
        function createFilterAndProject(obj)
            s.beta = 16.0;
            s.eta = 0.0;
            s.filterStep = 'LUMP';
            s.mesh = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.filter = FilterAndProject(s);
        end     

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            obj.computeDensity(x);  
            [f,dfdx]= obj.eigModes.computeFunctionAndGradient(obj.density);
            obj.value = f;  
            obj.gradient = dfdx;       
        end
                
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
        end
        
        function computeDensity(obj, x)  
            % Not Rounding Densities
%             s.operation = @(xV) obj.computeComplementaryDensity(x.fun,xV);
%             densComp = DomainFunction(s);
            
%             % Rounding Densities
%             s.operation = @(xV) obj.computeRoundedComplementaryDensity(x.fun,xV);
%             densComp = DomainFunction(s);

            % Filter and Project
            densComp = obj.filter.compute(1 - x.fun, 2);

             obj.density = densComp;
        end
        
        function rho = computeComplementaryDensity(obj,fun,xV)
            rho = fun.evaluate(xV);
            rho = 1 - rho;
        end

        function rho = computeRoundedComplementaryDensity(obj,fun,xV)
            rho = fun.evaluate(xV);
            rho = 1 - rho;
            rho = round(rho);
            rho = max(0,min(1,rho));
         end

    end
    
    
end