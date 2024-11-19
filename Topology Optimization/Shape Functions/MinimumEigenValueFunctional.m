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
       filterAdjoint
    end
    
    methods (Access = public)
        
        function obj = MinimumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            obj.computeDensity(x);  
%             obj.density.plot()
            [f,dfdx]= obj.eigModes.computeFunctionAndGradient(obj.density);    
            if ~isempty(obj.filterAdjoint)
                dfdx     = obj.filterAdjoint.compute(dfdx,2);
            end
            obj.value = f;  
            obj.gradient = dfdx;       
        end
                
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            if isfield(cParams,'filter')
                obj.filter         = cParams.filter;
            end
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint  = cParams.filterAdjoint;
            end
        end
        
        function x = filterDesignVariable(obj,x)
            x = obj.filter.compute(x,2);
            if ~isempty(obj.filterAdjoint) 
                obj.filterAdjoint.updateFilteredField(x);
            end
        end

        function computeDensity(obj,x)
            if isempty(obj.filter)
                densDomain  = x.fun;
                s.operation = @(xV) obj.computeComplementaryDensity(densDomain,xV);
                densHole = DomainFunction(s);
                obj.density = densHole;
            else
                xD  = 1 - x.obtainDomainFunction();
                xR = obj.filterDesignVariable(xD);
                obj.density = xR;
            end
         end
      
        function rho = computeComplementaryDensity(obj,fun,xV)
            rho = fun.evaluate(xV);
            rho = 1 - rho;
            rho = round(rho);
            rho = max(0,min(1,rho));
        end

    end
    
    
end