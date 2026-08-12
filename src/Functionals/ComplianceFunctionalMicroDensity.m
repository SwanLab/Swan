classdef ComplianceFunctionalMicroDensity < handle

    properties (Access = private)
        value0
        mesh
        filterB
        filterRho
        compliance
        material
        oldCost
        oldGradient
        xOld
    end

    methods (Access = public)

        function obj = ComplianceFunctionalMicroDensity(cParams)
            obj.init(cParams);
        end

        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            dx1 = xR{1} - obj.xOld{1};
            dx2 = xR{2} - obj.xOld{2};
            h = obj.mesh.computeMeanCellSize();
            V = obj.mesh.computeVolume();
            n1 = Norm(dx1,'H1',h)/V;
            n2 = Norm(dx2,'H1',h)/V;
            if max(n1,n2) > 0.005
                obj.material.setDesignVariable(xR);
                [J,dJ] = obj.computeComplianceFunctionAndGradient(x);
                obj.oldCost = J;
                obj.oldGradient = dJ;
                obj.xOld = xR;
            else
                sp1 = ScalarProduct(obj.oldGradient{1},dx1,'L2');
                sp2 = ScalarProduct(obj.oldGradient{2},dx2,'L2');
                J = obj.oldCost + sp1 + sp2;
                dJ = obj.oldGradient;
            end
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
            obj.xOld{1} = 1000;
            obj.xOld{2} = 1000;
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