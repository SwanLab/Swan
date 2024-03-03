classdef ComplianceFunctionalFromVademecum < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        compliance
        material
    end

    methods (Access = public)
        function obj = ComplianceFunctionalFromVademecum(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xR{1} = obj.filterDesignVariable(x.fun{1});
            xR{2} = obj.filterDesignVariable(x.fun{2});
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient();
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.material             = cParams.material;            
            obj.compliance = cParams.complainceFromConstitutive;
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,'QUADRATIC');
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();            
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            dJ     = obj.filter.compute(dJ,'QUADRATIC');
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
end