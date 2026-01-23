classdef ComplianceFunctionalThermoElastic < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        compliance
        material
        conductivity
        chiB
        kappaB
    end

    methods (Access = public)
        function obj = ComplianceFunctionalThermoElastic(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient(x,xR);
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.material   = cParams.material;
            obj.chiB       = cParams.chiB;
            obj.kappaB     = cParams.kappaB;
            obj.conductivity = cParams.conductivity;
            obj.compliance = cParams.complianceFromConstitutive;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj,x,xR)
            kappa  = obj.createDomainFunction(obj.conductivity.fun,xR);           % conductivity on the new domain
            dkappa = obj.createDomainFunction(obj.conductivity.dfun,xR); 
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            dC  = ChainRule.compute(x,dC);
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC,kappa,dkappa);
            dJ     = obj.filterFields(dJ);
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end

        function f = createDomainFunction(obj,fun,xR)
            s.operation = @(xV) obj.createConductivityAsDomainFunction(fun,xR{1},xV);
            s.mesh      = obj.mesh;
            kappa = DomainFunction(s);
            s.operation = @(xV) kappa.evaluate(xV).*(1-obj.chiB.evaluate(xV)) + obj.kappaB.*obj.chiB.evaluate(xV);
            s.mesh      = obj.mesh;
            f = DomainFunction(s);
        end

        function fV = createConductivityAsDomainFunction(obj,fun,xR,xV)
            densV = xR.evaluate(xV);
            fV = fun(densV);
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
            title = 'Compliance';
        end
    end
end