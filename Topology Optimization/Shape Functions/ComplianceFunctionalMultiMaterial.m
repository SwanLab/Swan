classdef ComplianceFunctionalMultiMaterial < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        compliance
        material
    end

    methods (Access = public)
        function obj = ComplianceFunctionalMultiMaterial(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            obj.material.setDesignVariable(xD);
            [J,dJ] = obj.computeComplianceFunctionAndGradient();
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.material   = cParams.material;
            obj.compliance = cParams.complainceFromConstitutive;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj)
            % AQUÍ ÉS ON ES FARIA EL DOBLE FOR I D'ON TREURIEM EL DJ TOTAL
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
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
            title = 'Compliance multimat';
        end
    end
end