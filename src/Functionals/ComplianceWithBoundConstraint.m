classdef ComplianceWithBoundConstraint < handle

    properties (Access = private)
        value0Compliance
    end

    properties (Access = private)
        mesh
        filterDesignVariable
        filterGradient
        compliance
        material
    end

    methods (Access = public)
        function obj = ComplianceWithBoundConstraint(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.density.obtainDomainFunction();
            xR = obj.filterDesignVariable.compute(xD{1},2);
            obj.filterGradient.updateFilteredField(xR);
            obj.material.setDesignVariable({xR});
            [Jc,dJc]   = obj.computeComplianceFunctionAndGradient();
            if isempty(obj.value0Compliance)
                obj.value0Compliance = Jc;
            end
            J          = Jc/obj.value0Compliance - x.bound;
            dJ.fValues = [dJc.fValues/obj.value0Compliance;-1];
            dJ         = {dJ};
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh                 = cParams.mesh;
            obj.filterDesignVariable = cParams.filterDesignVariable;
            obj.filterGradient       = cParams.filterGradient;
            obj.material             = cParams.material;
            obj.compliance           = cParams.complainceFromConstitutive;
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj)
            C      = obj.material.obtainTensor();
            dC     = obj.material.obtainTensorDerivative();
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            dJ     = obj.filterGradient.compute(dJ{1},2);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'ComplianceBoundConstr';
        end
    end
end