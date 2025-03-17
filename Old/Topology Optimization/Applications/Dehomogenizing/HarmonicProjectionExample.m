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
        unitBallProjector
    end

    methods (Access = public)

        function obj = HarmonicProjectionExample()
            obj.init();
            obj.loadDataExperiment();
            obj.createMesh();
            obj.createBoundaryMesh();
            obj.createHarmonicProjection();
            obj.createUnitBallProjector();
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


            theta(:,1) = atan2(alpha(:,1),alpha(:,2));  

            obj.plotOrientation(theta,1);
            alpha = obj.projectInUnitBall(alpha);
            theta(:,1) = atan2(alpha(:,1),alpha(:,2));
            obj.plotOrientation(theta,1);
            obj.orientationAngle = theta;
        end

        function vI = interpolateOrientationAngle(obj,v0)
            s.mesh    = obj.mesh;
            s.fValues = v0;
            s.order   = 'P0';
            p = LagrangianFunction(s);
            vI = p.projectToLinearNodalFunction();
        end

        function plotOrientation(obj,t,iFigure)
            figure(1)
            subplot(1,2,iFigure)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = cos(t);
            ty = sin(t);
            quiver(x,y,tx,ty)
            figure(2)
            subplot(1,2,iFigure)
            s.mesh = obj.mesh;
            s.field = t;
            p = NodalFieldPlotter(s);
            p.plot()
        end

        function project(obj)
            theta0 = obj.orientationAngle;
            u      = theta0;
            lambda = obj.computeInitialLambda();
            error = 1;
            i = 1;
            while error > 1e-12
                cost(i)      = obj.computeCost(u,theta0);
                optPrimal(i) = obj.computePrimalOptimaility(lambda,u,theta0);
                optDual(i)   = obj.computeDualHarmonicOptimality(u);
                error = norm([optPrimal,optDual]);
                [u,lambda]   = obj.solveProblem(u);

                figure(100)
                clf
                plot(cost,'-+')

                figure(101)
                clf
                hold on
                plot(optPrimal','-+')
                plot(optDual','-+')
                hold off

                obj.plotOrientation(u,2)
            end
        end

        function [v,lambda] = solveProblem(obj,vH)
            h  = obj.harmonicProjector;
            [v,lambda] = h.solveProblem(vH);
        end

        function createHarmonicProjection(obj)
            s.mesh = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;
            h = HarmonicProjector(s);
            obj.harmonicProjector = h;
        end

        function createUnitBallProjector(obj)
            u = UnitBallProjector([]);
            obj.unitBallProjector = u;
        end

        function vP = projectInUnitBall(obj,v)
            u = obj.unitBallProjector;
            vP = u.project(v);
        end

        function alpha = projectInHarmonicSpace(obj,alpha)
            h  = obj.harmonicProjector;
            nDim    = size(alpha,2);
            for iDim = 1:nDim
                vI  = alpha(:,iDim);
                vPI = h.project(vI);
                alpha(:,iDim) = vPI;
            end
        end

        function lambda0 = computeInitialLambda(obj)
            h = obj.harmonicProjector;
            lambda0 = h.computeInitalLambda();
        end

        function c = computeCost(obj,v,vH)
            h = obj.harmonicProjector;
            c = h.computeCost(v,vH);
        end

        function d = computePrimalOptimaility(obj,lambda,v,vH)
            h = obj.harmonicProjector;
            d = h.computePrimalOptimaility(lambda,v,vH);
        end

        function d = computeDualHarmonicOptimality(obj,v)
            h = obj.harmonicProjector;
            d = h.computeDualOptimality(v);
        end

        function dualOptT = computeDualUnitBallOptimality(obj,v)
            u = obj.unitBallProjector;
            dualOptT = u.computeDualOptimality(v);
        end

    end

end
