classdef PerimeterFunctional < handle

    properties (Access = private)
        mesh
        base
        domainFilter
        epsilon
        value0
    end

    properties (Access = private)
        riszFilter
        baseFun
    end

    methods (Access = public)
        function obj = PerimeterFunctional(cParams)
            obj.init(cParams);
            obj.domainFilter.updateEpsilon(obj.epsilon);
            obj.createRiszFilter();
            obj.createBaseFunction();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            J  = obj.computeFunction(xD{1},xR{1});
            dJ{1} = obj.computeGradient(xR{1});
            J  = obj.computeNonDimensionalValue(J);
            dJVal = obj.computeNonDimensionalValue(dJ{1}.fValues);
            dJ{1}.setFValues(dJVal);
        end

        function updateEpsilon(obj,epsilon)
            obj.domainFilter.updateEpsilon(epsilon);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.base         = cParams.uMesh;
            obj.domainFilter = cParams.filter;
            obj.epsilon      = cParams.epsilon;
            obj.value0       = cParams.value0;
        end

        function createRiszFilter(obj)
            s.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh         = obj.mesh;
            obj.riszFilter = FilterLump(s);
        end

        function createBaseFunction(obj)
            f           = CharacteristicFunction.create(obj.base);
            obj.baseFun = obj.riszFilter.compute(f,2);
        end

        function xR = filterDesignVariable(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.domainFilter.compute(x{i},2);
            end
        end

        function J = computeFunction(obj,xD,xR)
            b   = obj.baseFun;
            f   = xD.*(1-xR).*b;
            int = Integrator.compute(f,obj.mesh,2);
            J   = 2/(obj.epsilon)*int;
        end

        function dJH1 = computeGradient(obj,xR)
            b    = obj.baseFun;
            dJL2 = 2/(obj.epsilon)*(1-2*xR).*b;
            dJH1 = obj.riszFilter.compute(dJL2,2);
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