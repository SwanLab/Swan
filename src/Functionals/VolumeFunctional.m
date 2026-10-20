classdef VolumeFunctional < handle

    properties (Access = private)
        mesh
        base
        test
        volumeFieldIndex
    end

    properties (Access = private)
        baseFun
        totalVolume
    end

    methods (Access = public)

        function obj = VolumeFunctional(cParams)
            obj.init(cParams);
            obj.createBaseFunction();
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD      = x.obtainDomainFunction();
            nDesVar = numel(xD);
            iVolume = obj.volumeFieldIndex;

            if iVolume < 1 || iVolume > nDesVar
                error('VolumeFunctional:InvalidVolumeFieldIndex', ...
                    'volumeFieldIndex must be between 1 and %d.',nDesVar);
            end

            J  = obj.computeFunction(xD{iVolume});
            dJ = cell(nDesVar,1);

            for i = 1:nDesVar
                if i == iVolume
                    dJ{i} = obj.computeGradient();
                else
                    dJ{i} = obj.computeZeroGradient();
                end
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.base = cParams.uMesh;
            obj.test = cParams.test;

            if isfield(cParams,'volumeFieldIndex')
                obj.volumeFieldIndex = cParams.volumeFieldIndex;
            else
                obj.volumeFieldIndex = 1;
            end
        end

        function createBaseFunction(obj)
            s.trial     = LagrangianFunction.create(obj.mesh,1,obj.test.order);
            s.mesh      = obj.mesh;
            riszFilter  = FilterLump(s);
            f           = CharacteristicFunction.create(obj.base);
            obj.baseFun = riszFilter.compute(f,2);
        end

        function createTotalVolume(obj)
            V = Integrator.compute(obj.baseFun,obj.mesh,2);
            obj.totalVolume = V;
        end

        function J = computeFunction(obj,x)
            volume = Integrator.compute(x.*obj.baseFun,obj.mesh,2);
            J      = volume/obj.totalVolume;
        end

        function dJ = computeGradient(obj)
            dJ = copy(obj.baseFun);
            dJ.setFValues(dJ.fValues./obj.totalVolume);
        end

        function dJ = computeZeroGradient(obj)
            dJ = copy(obj.baseFun);
            dJ.setFValues(zeros(size(dJ.fValues)));
        end

    end

    methods (Static, Access = public)

        function title = getTitleToPlot()
            title = 'Volume';
        end

    end

end