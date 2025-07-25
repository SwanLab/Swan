classdef ThermalComplianceFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        compliance
        material
        stateProblem
    end

    methods (Access = public)
        function obj = ThermalComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            % make sure that the conductivity is being updated!---
            u  = obj.computeStateVariable();
            J = obj.computeFunction(u);
            dJ = obj.computeGradient(u);
        end

    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.stateProblem = cParams.stateProblem;
%             obj.material   = cParams.material; ----
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

        function u = computeStateVariable(obj)
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
        end

        function J = computeFunction(obj, u)
            s.operation = @(xV) evaluate(u,xV);
            s.mesh = u.mesh;
            dom = DomainFunction(s);
            dCompliance = dom;
            J  = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);    
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
        end

        function fVR = evaluate(u,xR,xV)
            kappa = obj.createDomainFunction(obj.conductivity.fun, xR); 
            thE = DP(kappa*Grad(u),Grad(u)); % check dimensions
            fVR = thE.evaluate(xV);
        end

        function dJ = computeGradient(obj,xR,u)
            dkappa = obj.createDomainFunction(obj.conductivity.dfun,xR); 
            dJ = - DP(dkappa*Grad(u),Grad(u));
            dJ     = obj.filterFields(dJ);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end

        function f = createDomainFunction(obj,fun,xR)
            s.operation = @(xV) obj.createConductivityAsDomainFunction(fun,xR,xV);
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