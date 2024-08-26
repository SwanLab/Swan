classdef PerimeterFunctional < handle

    properties (Access = private)
        mesh
        filter
        epsilon
        value0
    end

    methods (Access = public)
        function obj = PerimeterFunctional(cParams)
            obj.init(cParams);
            obj.filter.updateEpsilon(obj.epsilon);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            J  = obj.computeFunction(xD,xR);
            dJ = obj.computeGradient(xR);
            J  = obj.computeNonDimensionalValue(J);
            dJ.fValues = obj.computeNonDimensionalValue(dJ.fValues);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.filter  = cParams.filter;
            obj.epsilon = cParams.epsilon;
            obj.value0  = cParams.value0;
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
        end

        function J = computeFunction(obj,xD,xR)
            f   = xD.*(1-xR);
            int = Integrator.compute(f,obj.mesh,2);
            J   = 2/(obj.epsilon)*int;
        end

        function dJ = computeGradient(obj,xR)
            dj        = 2/(obj.epsilon)*(1-2*xR.fValues);
            s.fValues = dj;
            s.mesh    = xR.mesh;
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