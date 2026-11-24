classdef LatticeVolumeConstraint < handle

    properties (Access = private)
        mesh
        volumeTarget
        sideFields
        totalVolume
    end

    methods (Access = public)

        function obj = LatticeVolumeConstraint(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            m1 = xD{obj.sideFields(1)};
            m2 = xD{obj.sideFields(2)};
            V  = obj.computeVolumeFraction(m1,m2);
            J  = V/obj.volumeTarget - 1;
            dJ = obj.computeGradient(xD,m1,m2);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.sideFields   = cParams.sideFields;
            obj.totalVolume  = obj.mesh.computeVolume();
        end

        function V = computeVolumeFraction(obj,m1,m2)
            hole = m1.*m2;
            V    = 1 - Integrator.compute(hole,obj.mesh,2)/obj.totalVolume;
        end

        function dJ = computeGradient(obj,xD,m1,m2)
            scale = obj.volumeTarget*obj.totalVolume;
            dJ    = cell(numel(xD),1);
            for i = 1:numel(xD)
                dJ{i} = copy(xD{i});
                dJ{i}.setFValues(zeros(size(xD{i}.fValues)));
            end
            dJ{obj.sideFields(1)}.setFValues(-m2.fValues/scale);
            dJ{obj.sideFields(2)}.setFValues(-m1.fValues/scale);
        end

    end

    methods (Static, Access = public)

        function title = getTitleToPlot()
            title = 'Volume constraint';
        end

    end

end