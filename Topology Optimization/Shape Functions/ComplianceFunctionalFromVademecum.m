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
            mat  = obj.computeMaterial(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient(mat);
        end
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.compliance = cParams.complainceFromConstitutive;
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,'LINEAR');
        end

        function C = computeMaterial(obj,x)
            s.type         = 'HomogenizedMicrostructure';
            s.fileName     = 'Rectangle';
            s.microParams = x;
            C  = Material.create(s);            
        end

        function dC = computeMaterialDerivative(obj,x)
            mI = obj.material;
            dC = mI.computeConstitutiveTensorDerivative(x);
        end        

        function [J,dJ] = computeComplianceFunctionAndGradient(obj,mat)
            [J,dJ] = obj.compliance.computeFunctionAndGradient(mat);
            dJ     = obj.filter.compute(dJ,'LINEAR');
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