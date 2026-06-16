classdef ComplianceFunctionalMicroDensity < handle

    properties (Access = private)
        value0
        mesh
        filterB
        filterRho
        compliance
        material
    end

    methods (Access = public)

        function obj = ComplianceFunctionalMicroDensity(cParams)
            obj.init(cParams);
        end

        function [J, dJ] = computeFunctionAndGradient(obj, x)
            xD = x.obtainDomainFunction();    
            xR = obj.filterFields(xD);        
            obj.material.setDesignVariable(xR);
            [J, dJ] = obj.computeComplianceFunctionAndGradient(x);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh       = cParams.mesh;
            obj.filterB    = cParams.filter;         
            obj.filterRho  = cParams.filterDensity;  
            obj.material   = cParams.material;
            obj.compliance = cParams.complainceFromConstitutive;
            if isfield(cParams, 'value0')
                obj.value0 = cParams.value0;
            end
        end

        function xR = filterFields(obj, x)
            % x{1} = b,   
            xR{1} = obj.filterB.compute(x{1}, 2);
            xR{2} = obj.filterRho.compute(x{2}, 2);
        end

        function [J, dJ] = computeComplianceFunctionAndGradient(obj, x)
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();  
            dC = ChainRule.compute(x, dC);
            [J, dJ] = obj.compliance.computeFunctionAndGradient(C, dC);
            dJ = obj.filterFields(dJ);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end

        function x = computeNonDimensionalValue(obj, x)
            x = x / obj.value0;
        end

        function dx = computeNonDimensionalGradient(obj, dx)
            for i = 1:length(dx)
                dx{i}.setFValues(dx{i}.fValues / obj.value0);
            end
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end

end