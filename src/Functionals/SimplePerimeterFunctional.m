classdef SimplePerimeterFunctional < handle

    properties (Access = private)
        mesh
        filter
        quadrature
        weight
        updated
        iter
        value0
    end

    methods (Access = public)
        function obj = SimplePerimeterFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
%             iter = x{2};
%             x = x{1};
            xD = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            J  = obj.computeFunction(xD{1},xR{1});
            dJ{1} = obj.computeGradient(xR{1});
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.updated = 0;
            obj.mesh    = cParams.mesh;
            obj.filter  = cParams.filter;
            obj.iter = 0;
        end

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end

        function J = computeFunction(obj, rho, rhoF)
            perimeter = rho.*(1-rhoF);
            J      = Integrator.compute(perimeter,obj.mesh,obj.quadrature.order);
        end

        function dJ = computeGradient(obj, rhoF)
            dj        = 1-2*rhoF.fValues;
            s.fValues = dj;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            dJ        = LagrangianFunction(s);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter';
        end
    end
end