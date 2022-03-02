classdef HarmonicProjectionExample < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        filePath
        fileName
        iteration
        experimentData
        mesh
        boundaryMesh
        orientationAngle
        orientationAngleGauss
        harmonicProjector
    end

    methods (Access = public)

        function obj = HarmonicProjectionExample()
            obj.init();
            obj.loadDataExperiment();
            obj.createMesh();
            obj.createBoundaryMesh();
            obj.createHarmonicProjection();
            obj.storeOrientationAngle();
            obj.project()
        end

    end

    methods (Access = private)

        function init(obj)
            obj.filePath = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
            obj.fileName = 'LshapeCoarseSuperEllipseDesignVariable';
            obj.iteration = 665;
        end

        function loadDataExperiment(obj)
            s.fileName = [obj.fileName,num2str(obj.iteration)];
            s.folderPath = fullfile(obj.filePath );
            w = WrapperMshResFiles(s);
            w.compute();
            obj.experimentData = w;
        end

        function createMesh(obj)
            d = obj.experimentData;
            obj.mesh = d.mesh;
        end

        function createBoundaryMesh(obj)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            b = boundary(x,y,1);
            obj.boundaryMesh = b;
        end

        function storeOrientationAngle(obj)
            d = obj.experimentData;
            alpha0  = d.dataRes.AlphaGauss;
            alpha(:,1) = obj.interpolateOrientationAngle(alpha0(:,1));
            alpha(:,2) = obj.interpolateOrientationAngle(alpha0(:,2));
            alpha = obj.normalizeVector(alpha);
            obj.orientationAngle = alpha;
        end

        function vI = interpolateOrientationAngle(obj,v0)
            s.mesh    = obj.mesh;
            s.fValues = v0;
            p = PieceWiseConstantFunction(s);
            vI = p.projectToLinearNodalFunction(); 
        end

        function vN = normalizeVector(obj,v)
            vx = v(:,1);
            vy = v(:,2);
            norm = sqrt(vx.^2 + vy.^2);
            vN(:,1) = vx./norm;
            vN(:,2) = vy./norm;
        end

        function plotOrientation(obj,t)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = t(:,1);
            ty = t(:,2);
            quiver(x,y,tx,ty)
        end

        function project(obj)
            alpha = obj.orientationAngle;
            obj.plotOrientation(alpha)
            w = zeros(size(alpha));
            z = alpha;
            for i = 1:100
                [u,errF1,errF2] = obj.harmonicProjection(z+w);
                z = obj.normalizeVector(u-w);
                w = w + z -u;
                if mod(i,1) == 0
                    obj.plotOrientation(z)
                    %
                end
            end
            alpha = obj.normalizeVector(alpha);
            obj.plotOrientation(alpha)
            plot(1:i,log([errF1' errF2']))
        end

        function createHarmonicProjection(obj)
            s.mesh = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;
            h = HarmonicProjection(s);
            obj.harmonicProjector = h;
        end

        function [alpha,errF1,errF2] = harmonicProjection(obj,alpha)
            h = obj.harmonicProjector;
            [alpha(:,1),errF1] = h.project(alpha(:,1));
            [alpha(:,2),errF2] = h.project(alpha(:,2));
        end

    end

end