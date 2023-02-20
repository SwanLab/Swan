classdef SingularitiesFinder < handle

    properties (Access = public)

    end

    properties (Access = private)
        orientation
        mesh
    end

    properties (Access = private)
        meshDisc
        isElemSingular
    end

    methods (Access = public)

        function obj = SingularitiesFinder(cParams)
            obj.init(cParams)
        end

        function isS = computeSingularElements(obj)
            obj.computeSingularities();
            isS = obj.isElemSingular;
        end

        function plot(obj)
            obj.plotOrientationVector();
            obj.plotSingularities();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation = cParams.orientation;
        end

        function computeSingularities(obj)
            a = obj.orientation;
            a = obj.mapP1ToP1Discontinous(a);
            a1 = a(:,:,1);
            a2 = a(:,:,2);
            a3 = a(:,:,3);
            a1a2 = obj.scalarProduct(a1,a2);
            a1a3 = obj.scalarProduct(a1,a3);
            a2a3 = obj.scalarProduct(a2,a3);
            isS = sign(a1a2.*a1a3.*a2a3);
            obj.isElemSingular = isS<0;
        end

        function plotOrientationVector(obj)
            figure()
            a = obj.orientation;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            ax = a(:,1);
            ay = a(:,2);
            q = quiver(x,y,ax,ay);
            q.ShowArrowHead = 'off';
        end

        function plotSingularities(obj)
            s.fValues = double(obj.isElemSingular);
            s.mesh    = obj.mesh;
            p0 = P0Function(s);
            p0.plot();
        end

        function fP1 = mapP1ToP1Discontinous(obj,f)
            s.fValues = f;
            s.mesh    = obj.mesh;
            fun = P1Function(s);
            funP1D = fun.project('P1D');
            fP1 = permute(funP1D.fValues, [3 1 2]);
        end

    end

    methods (Access = private, Static)

        function ab = scalarProduct(a,b)
            ab = a(:,1).*b(:,1) + a(:,2).*b(:,2);
        end

    end

end