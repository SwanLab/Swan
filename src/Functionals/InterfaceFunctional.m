classdef InterfaceFunctional < handle

    properties (Access = private)
        mesh
        base
        filter
        epsilon
        value0
        sign
        signF
        ds
        incSign
        valueOld
        tarVolume
    end

    properties (Access = private)
        riszFilter
        baseFun
        totalVolume
    end

    methods (Access = public)
        function obj = InterfaceFunctional(cParams)
            obj.init(cParams);
            obj.filter.updateEpsilon(obj.epsilon);
            obj.createRiszFilter();
            obj.createBaseFunction();
            obj.computeTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            V  = Integrator.compute(xD{1},x.fun.mesh,2);
            V  = V/(obj.tarVolume*obj.totalVolume) - 1;
            obj.updateSignParameter(xD{1},V);
            xR    = obj.filterDesignVariable(xD);
            J     = obj.computeFunction(xD{1},xR{1});
            dJ{1} = obj.computeGradient(xR{1});
            J     = obj.computeNonDimensionalValue(J);
            dJVal = obj.computeNonDimensionalValue(dJ{1}.fValues);
            dJ{1}.setFValues(dJVal);
        end

        function updateEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
            obj.filter.updateEpsilon(epsilon);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.base      = cParams.uMesh;
            obj.filter    = cParams.filter;
            obj.epsilon   = cParams.epsilon;
            obj.value0    = cParams.value0;
            obj.sign      = cParams.signInitial;
            obj.signF     = cParams.signFinal;
            obj.ds        = (obj.signF-obj.sign)/40;
            obj.incSign   = false;
            obj.valueOld  = -inf;
            obj.tarVolume = cParams.tarVolume;
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

        function computeTotalVolume(obj)
            u = ConstantFunction.create(1,obj.mesh);
            obj.totalVolume = Integrator.compute(u,obj.mesh,2);
        end

        function xR = filterDesignVariable(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function J = computeFunction(obj,xD,xR)
            b   = obj.baseFun;
            f   = xD.*(1-xR).*b;
            int = Integrator.compute(f,obj.mesh,2);
            J   = 2*obj.sign/(obj.epsilon)*int;
        end

        function dJH1 = computeGradient(obj,xR)
            b    = obj.baseFun;
            dJL2 = 2*obj.sign/(obj.epsilon)*(1-2*xR).*b;
            dJH1 = obj.riszFilter.compute(dJL2,2);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end

        function updateSignParameter(obj,x,V)
            delta = norm(x.fValues-obj.valueOld);
            if ~obj.incSign
                if V<0.1
                    obj.incSign = true;
                end
            elseif obj.incSign && delta==0
                obj.sign = min(obj.sign+obj.ds,obj.signF);
            end
            obj.valueOld = x.fValues;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter';
        end
    end
end