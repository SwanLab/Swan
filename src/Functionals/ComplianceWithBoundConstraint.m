classdef ComplianceWithBoundConstraint < handle

    properties (Access = private)
        value0Compliance
    end

    properties (Access = private)
        mesh
        filterDesignVariable
        filterGradient
        compliance
        C
        dC
    end

    methods (Access = public)
        function obj = ComplianceWithBoundConstraint(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.density.obtainDomainFunction();
            xR = obj.filterDesignVariable.compute(xD{1},2);
            obj.filterGradient.updateFilteredField(xR);
            [Jc,dJc]   = obj.computeComplianceFunctionAndGradient(xR);
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
            obj.C                    = cParams.C;
            obj.dC                   = cParams.dC;
            obj.compliance           = cParams.complainceFromConstitutive;
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj,xR)
            CxR     = obj.C(xR);
            dCxR{1} = obj.dC(xR);
            [J,dJ]  = obj.compliance.computeFunctionAndGradient(CxR,dCxR);
            dJ      = obj.filterGradient.compute(dJ{1},2);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'ComplianceBoundConstr';
        end
    end
end