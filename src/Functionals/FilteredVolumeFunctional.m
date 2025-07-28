classdef FilteredVolumeFunctional < handle

    properties (Access = private)
        mesh
        base
        filter
    end

    properties (Access = private)
        baseFun
        totalVolume
    end

    methods (Access = public)
        function obj = FilteredVolumeFunctional(cParams)
            obj.init(cParams);
            obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            J  = obj.computeFunction(xR{1});
            dJ{1} = obj.computeGradient();
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.base   = cParams.uMesh;
            obj.filter = cParams.filter;
        end

        function createBaseFunction(obj)
            obj.baseFun = CharacteristicFunction.create(obj.base);
        end

        function createTotalVolume(obj)
            dV = obj.baseFun;
            V  = Integrator.compute(dV,obj.mesh,2);
            obj.totalVolume = V;
        end

        function xR = filterDesignVariable(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end

        function J = computeFunction(obj,xR)
            b   = obj.baseFun;
            f   = b.*xR./obj.totalVolume;
            int = Integrator.compute(f,obj.mesh,2);
            J   = int;
        end

        function dJ = computeGradient(obj)
            b  = obj.baseFun;
            dj = obj.filter.compute(b,2);
            dJ = dj./obj.totalVolume;
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Filtered volume';
        end
    end
end