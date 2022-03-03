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
            obj.plotOrientation(alpha,1);
            alpha = obj.projectInUnitBall(alpha);
            obj.plotOrientation(alpha,2);
            obj.orientationAngle = alpha;
        end

        function vI = interpolateOrientationAngle(obj,v0)
            s.mesh    = obj.mesh;
            s.fValues = v0;
            p = PieceWiseConstantFunction(s);
            vI = p.projectToLinearNodalFunction(); 
        end

        function plotOrientation(obj,t,iFigure)
            figure(1)
            subplot(1,3,iFigure)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = t(:,1);
            ty = t(:,2);
            quiver(x,y,tx,ty)
        end

        function project(obj)
            alpha0 = obj.orientationAngle;            
            z = zeros(size(alpha0));
            v = alpha0;
            t = 1*1e-2;
            for i = 1:100000
                j = 2*(i-1)+1;
                
                r  = (0*alpha0+t*(v-z))/(0+t);
                u = obj.projectInHarmonicSpace(r);
                cost(j)     = obj.computeCost(u,alpha0);
                optDualH(j) = obj.computeDualHarmonicOptimality(u);
                optDualU(j) = obj.computeDualUnitBallOptimality(u);

               
                v = obj.projectInUnitBall(u+z);
                cost(j+1)     = obj.computeCost(v,alpha0);
                optDualH(j+1) = obj.computeDualHarmonicOptimality(v);
                optDualU(j+1) = obj.computeDualUnitBallOptimality(v);                





%                 figure(102)
%                 clf
%                 hold on
%                 plot(optDualH(1:2:j),'-+')
%                 plot(optDualU(1:2:j),'-+')
%                 hold off

                if mod(i,100) == 0      
                figure(100)
                clf
                hold on
                plot(1:2:j,cost(1:2:j)','-+')                
                plot(2:2:j+1,cost(2:2:j+1)','-+')
                hold off

                figure(101)
                clf
                hold on
                plot(1:2:j,optDualU(1:2:j)','-+')                
                plot(2:2:j+1,optDualH(2:2:j+1)','-+')
                hold off
                 obj.plotOrientation(u,1) 
                 obj.plotOrientation(v,2)     
                 obj.plotOrientation(z,3)                  
                end
                z = z + t*(u-v);               

            end
            alpha0 = obj.projectInUnitBall(alpha0);
            obj.plotOrientation(alpha0)
            plot(log([cost' optDualH' optDualU']))
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

        function costT = computeCost(obj,v,vP)
            h     = obj.harmonicProjector;            
            nDim  = size(v,2);
            cost  = zeros(nDim,1);            
            for iDim = 1:nDim
                vI  = v(:,iDim);
                vPI = vP(:,iDim);
                c   = h.computeCost(vI,vPI);             
                cost(iDim) = c;
            end            
            costT = sum(cost);
        end        

        function dualOptT = computeDualHarmonicOptimality(obj,v)
            h       = obj.harmonicProjector;            
            nDim    = size(v,2);
            dualOpt = zeros(nDim,1);            
            for iDim = 1:nDim
                vI  = v(:,iDim);
                d   = h.computeDualOptimality(vI);             
                dualOpt(iDim) = d;
            end            
            dualOptT = sum(dualOpt);
        end   

        function dualOptT = computeDualUnitBallOptimality(obj,v)
            u = obj.unitBallProjector;
            dualOptT = u.computeDualOptimality(v);
        end             

    end

end