classdef MinimumEigenValueFunctional < handle
    
    properties (Access = public)
        value
    end
       
    properties (Access = private)
       eigModes 
       designVariable
       mesh
       filter
       filterAdjoint
       isComplementary
    end
    
    methods (Access = public)
        
        function obj = MinimumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            xD  = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD{1});
            if obj.isComplementary == true
                xR.setFValues(max(min(1 - xR.fValues,1),0)); % 1 - FP
            end
            [f,dfdx]= obj.eigModes.computeFunctionAndGradient(xR);    
            if ~isempty(obj.filterAdjoint)
                dfdx     = obj.filterAdjoint.compute(dfdx,2);
            else
                dfdx     = obj.filter.compute(dfdx,2);
            end
            if obj.isComplementary == true
                dfdx.setFValues(-dfdx.fValues);              % Chain rule for (1 - FP)
            end
        end
                
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            obj.isComplementary  = cParams.isComplementary;
            obj.filter = cParams.filter;
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint  = cParams.filterAdjoint;
            end
        end
      
        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
            if ~isempty(obj.filterAdjoint)
                obj.filterAdjoint.updateFilteredField(xR);
            end
        end

    end
    
    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'MinimumEigenvalue';
        end
    end  
end