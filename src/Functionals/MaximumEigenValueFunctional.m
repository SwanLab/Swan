classdef MaximumEigenValueFunctional < handle
       
    properties (Access = private)
       density
       eigModes 
       designVariable
       mesh
       filter
       filterAdjoint
       iter
       value0
    end
    
    methods (Access = public)
        
        function obj = MaximumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [lambdas, phis] = computeEigenModes(obj, x, n)
            xD  = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD{1});     
            xR.setFValues(1 - xR.fValues);
            [lambdas, phis] = obj.eigModes.getEigenModesComputer(xR,n);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = StiffnessEigenModesDisplacementComputer(cParams);
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            obj.filter = cParams.filter;
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint = cParams.filterAdjoint;
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
                xFiltered = obj.filter.getFilteredField();
                obj.filterAdjoint.updateFilteredField(xFiltered);
            end
        end
    end
end