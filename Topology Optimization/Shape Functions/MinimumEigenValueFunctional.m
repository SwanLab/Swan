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
        
        function computeFunctionAndGradient(obj)            
            obj.computeFunction();
            obj.computeGradient();
        end
        
        function computeFunction(obj)
            obj.computeDensity();
            obj.value = obj.eigModes.compute(obj.density);
            %obj.normalizeFunction();
        end
        
        function computeGradient(obj)
            obj.computeDensity();
            dfdx = obj.eigModes.computeGradient(eigN);
            obj.gradient = dfdx';
            %obj.normalizeGradient();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
        end
        
        function computeDensity(obj)
            dens  = obj.designVariable.fun.project('P0');            
            dens.fValues = 1 - dens.fValues; 
            dens.fValues = round(dens.fValues);   
            obj.density = dens;
        end
        
    end
    
    
end