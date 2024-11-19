classdef ComplianceFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        filterAdjoint
        compliance
        material
    end

    methods (Access = public)
        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient();
        end

        function updateFilterParams(obj, beta)
            obj.filter.updateBeta(beta);
            obj.filterAdjoint.updateBeta(beta);
        end
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.material   = cParams.material;
            obj.compliance = cParams.complainceFromConstitutive;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint = cParams.filterAdjoint;            
            end
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
            if ~isempty(obj.filterAdjoint) 
                obj.filterAdjoint.updateFilteredField(xR);
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            if isempty(obj.filterAdjoint) 
                dJ     = obj.filter.compute(dJ,2);
            else
                dJ     = obj.filterAdjoint.compute(dJ,2);
            end

            if isempty(obj.value0)
                obj.value0 = J;
            end
            J          = obj.computeNonDimensionalValue(J);
            dJ.fValues = obj.computeNonDimensionalValue(dJ.fValues);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end