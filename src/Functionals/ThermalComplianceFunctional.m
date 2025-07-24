classdef ThermalComplianceFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        filterAdjoint
        compliance
        material
        iter
        filteredDesignVariable
    end

    methods (Access = public)
        function obj = ThermalComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeThermalComplianceFunctionAndGradient(x);
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.material   = cParams.material;
        end

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [J,dJ] = computeThermalComplianceFunctionAndGradient(obj,x)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            dC  = ChainRule.compute(x,dC);
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end

        function J = computeFunction(obj,C,u)
            dCompliance = ElasticEnergyDensity(C,u);
            J           = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);
        end

        function dj = computeGradient(dkappa,u)
            nDesVar = length(dC);
            dj      = cell(nDesVar,1);
            for i = 1:nDesVar
                dj{i}   = - dkappa{i}.*DP(Grad(u), Grad(u));
            end
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end

        function dx = computeNonDimensionalGradient(obj,dx)
            refX = obj.value0;
            for i = 1:length(dx)
                dx{i}.setFValues(dx{i}.fValues/refX);
            end
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Thermal Compliance';
        end
    end
end