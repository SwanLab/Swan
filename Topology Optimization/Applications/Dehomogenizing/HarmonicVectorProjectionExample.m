classdef HarmonicVectorProjectionExample < handle

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
        orientationVector
    end

    methods (Access = public)

        function obj = HarmonicVectorProjectionExample()
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
          %  obj.plotAlpha0(alpha0);
            alpha(:,1) = obj.interpolateOrientationAngle(alpha0(:,1));
            alpha(:,2) = obj.interpolateOrientationAngle(alpha0(:,2));
           
            x2 = obj.mesh.coord(:,2);
            x1 = obj.mesh.coord(:,1);
            
            %theta = atan2(x2,x1);            
            %thetaV = pi/6;
%            theta = thetaV*ones(size(x2));
            
           % theta = rand(1)*x1 + rand(1)*x2 + rand(1);


          %  alpha(:,1) = cos(theta);
          %  alpha(:,2) = sin(theta);            
            
            
            theta(:,1) = atan2(alpha(:,1),alpha(:,2));  
            obj.plotOrientation(alpha,1);
            alpha = obj.projectInUnitBall(alpha);
            theta(:,1) = atan2(alpha(:,1),alpha(:,2));  
            obj.plotOrientation(alpha,1);
            obj.orientationAngle = theta;
            obj.orientationVector = alpha;
        end

        function plotAlpha0(obj,t)
            m    = obj.mesh;
            coord0 = m.computeBaricenter';
            x = coord0(:,1);
            y = coord0(:,2);
            tx = t(:,1);
            ty = t(:,2);     
            figure(5)
            quiver(x,y,tx,ty)
        end

        function vI = interpolateOrientationAngle(obj,v0)
            s.mesh    = obj.mesh;
            s.fValues = v0;
            p = PieceWiseConstantFunction(s);
            vI = p.projectToLinearNodalFunction(); 
        end

        function plotOrientation(obj,t,iFigure)
            figure(1)
            subplot(1,2,iFigure)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = t(:,1);
            ty = t(:,2);
            quiver(x,y,tx,ty)
            figure(2)
            subplot(1,2,iFigure)
            s.mesh = obj.mesh;
            s.field = atan2(t(:,1),t(:,2));
            p = NodalFieldPlotter(s);
            p.plot()
        end

        function project(obj)
            alpha0 = obj.orientationVector;            

            figure(23)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = alpha0(:,1);
            ty = alpha0(:,2);
            quiver(x,y,tx,ty)

            theta = atan2(alpha0(:,1),alpha0(:,2));
            dalpha0(:,1) = cos(2*theta);
            dalpha0(:,2) = sin(2*theta);


            figure(24)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = dalpha0(:,1);
            ty = dalpha0(:,2);
            quiver(x,y,tx,ty)


            u      = alpha0;
            lambda = obj.computeInitialLambda();
            isErrorLarge = true;
            i = 1;
            while isErrorLarge
                cost(i)      = obj.computeCost(u,alpha0);
              %  optPrimal(i) = obj.computePrimalOptimaility(lambda,u,alpha0);
              %  optDual(i)   = obj.computeDualHarmonicOptimality(u);
              %  error = norm([optPrimal,optDual]);
                
                theta = atan2(u(:,1),u(:,2));
                u2(:,1) = cos(2*theta);
                u2(:,2) = sin(2*theta);
                [uNew2,lambda] = obj.solveProblem(dalpha0,u2);
                theta = 0.5*atan2(uNew2(:,1),uNew2(:,2));
                uNew(:,1) = cos(theta);
                uNew(:,2) = sin(theta);
                
                
                
                uNew = obj.projectInUnitBall(uNew);
                err(i) = norm(uNew(:)-u(:))/norm(u(:));
                u = uNew;
                figure(200)
                clf
                plot(cost,'-+')                

                figure(201)
                plot(err,'-+')    
% 
%                 figure(101)
%                 clf
%                 hold on
%                 plot(optPrimal','-+')                
%                 plot(optDual','-+')
%                 hold off
                isErrorLarge = err(i) > 1e-3;
                i = i +1;
                obj.plotOrientation(u,2) 
       

            end
        end

        function [v,lambda] = solveProblem(obj,alpha0,vH)
           h  = obj.harmonicProjector;
           [v,lambda] = h.solveProblem(alpha0,vH);
        end

        function createHarmonicProjection(obj)
            s.mesh = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;
            h = LinearizedHarmonicProjector(s);
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