classdef FilteredVolumeFunctional < handle

    properties (Access = private)
        mesh
        base
        filter
        p
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
            Vp = obj.computeVolume(xR{1});
            J  = obj.computeFunction(Vp);
            dJ{1} = obj.computeGradient(xR{1},Vp);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.base   = cParams.uMesh;
            obj.filter = cParams.filter;
            obj.p      = cParams.p;
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

        function Vp = computeVolume(obj,xR)
            b   = obj.baseFun;
            xP  = xR.^obj.p;
            int = Integrator.compute(b.*xP,obj.mesh,3);
            Vp  = int^(1/obj.p);
        end

        function J = computeFunction(obj,Vp)
            J = Vp/(obj.totalVolume^(1/obj.p));
        end

        function dJ = computeGradient(obj,xR,Vp)
            b  = obj.baseFun;
            dj = ((Vp^(1-obj.p))*xR.^(obj.p-1))./(obj.totalVolume^(1/obj.p));
            dJ = obj.filter.compute(b.*dj,3);
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Filtered volume';
        end
    end
end