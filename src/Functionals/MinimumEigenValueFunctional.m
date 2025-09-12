classdef MinimumEigenValueFunctional < handle
    
    properties (Access = public)
        value
    end
       
    properties (Access = private)
       density
       eigModes 
       designVariable
       mesh
       filter
       filterAdjoint
       iter
       isCompl
    end
    
    methods (Access = public)
        
        function obj = MinimumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            xD  = x.obtainDomainFunction();             % rho
            xR = obj.filterDesignVariable(xD{1});       % FP rho
            if obj.isCompl == true
                xR.setFValues(1 - xR.fValues);          % 1 - FP
            end
            [f,dfdx]= obj.eigModes.computeFunctionAndGradient(xR);    
            if ~isempty(obj.filterAdjoint)
                dfdx     = obj.filterAdjoint.compute(dfdx,2);
            else
                if isa(obj.filter, 'HeavisideProjector')
                    sensitVals         = obj.filter.derive(xR);
                    dfdx     = dfdx.project('P1',obj.mesh);
                    dfdx.fValues       = dfdx.fValues.*sensitVals;
                else
                    dfdx     = obj.filter.compute(dfdx,2);
                end
            end
            if obj.isCompl == true
                dfdx.setFValues(-dfdx.fValues);          % Chain rule
            end
        end
                
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            if isfield(cParams,'filter')
                obj.filter = cParams.filter;
            end
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint  = cParams.filterAdjoint;
            end
            if isfield(cParams,'isCompl')
                obj.isCompl  = cParams.isCompl;
            end
        end
      
        function xR = filterDesignVariable(obj,x)
            if isa(obj.filter, 'HeavisideProjector')
                fValues = obj.filter.project(x);
                xR = FeFunction.create(x.order,fValues,obj.mesh);
            else
                xR = obj.filter.compute(x,2);
            end
            if ~isempty(obj.filterAdjoint)
                obj.filterAdjoint.updateFilteredField(obj.filter);
            end
        end

    end
    
    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'MinimumEigenvalue';
        end
    end  
end