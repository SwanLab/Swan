classdef PerimeterFunctionalParameter < handle

    properties (Access = private)
        mesh
        perimeter
        dJ     
    end

    properties (Access = private)
    end

    methods (Access = public)
        function obj = PerimeterFunctionalParameter(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD    = x.obtainDomainFunction();
            J     = obj.computeFunction(xD{1});
            dJ{1} = obj.computeGradient(xD{1});
            J     = obj.computeNonDimensionalValue(J);
            dJVal = obj.computeNonDimensionalValue(dJ{1}.fValues);
            dJ{1}.setFValues(dJVal);
        end
   

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.perimeter = cParams.perimeter;
            obj.dJ        = cParams.dJ;
        end

        function J = computeFunction(obj,xD)
            J = obj.perimeter(xD);
        end

        function dJ = computeGradient(obj,xD)
            dJ = obj.dJ(xD);
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