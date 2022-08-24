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
        backgroundMesh
        superEllipse
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
        %    obj.dehomogenize();            
            obj.project()
            obj.dehomogenize();                        
        end

    end

    methods (Access = private)

        function init(obj)
            close all
            obj.filePath = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
            obj.fileName = 'LshapeCoarseSuperEllipseDesignVariable';
            obj.iteration = 665;

%             obj.filePath = '/home/alex/git-repos/Swan/Topology Optimization/Applications/Dehomogenizing/ExampleCompliance/';  
%             obj.fileName = 'ExperimentingPlotSuperEllipse';
%             obj.iteration = 64;
                        
        end

        function loadDataExperiment(obj)
            s.fileName = [obj.fileName,num2str(obj.iteration)];
            s.folderPath = fullfile(obj.filePath );
            w = WrapperMshResFiles(s);
            w.compute();
            obj.experimentData = w;
        end

        function createMesh(obj)
            %d = obj.experimentData;
            %obj.mesh = d.mesh;          
         %   obj.mesh = obj.createSquareMesh();
            obj.mesh  = obj.createCircularMesh();
        end

        function m = createSquareMesh(obj)
            h = 0.05;
            [X,Y] = meshgrid(-1:h:1,0:h:1); 
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
            s.connec = delaunay(s.coord); 
            m = Mesh(s);            
            m.plot;
        end

        function m = createCircularMesh(obj)
             %h = 0.03;
%             [X,Y] = meshgrid(-1:h:1,0:h:1);                        
%             isD = (X(:)-0).^2 + (Y(:)-0).^2 - 0.15^2 > 0;
%             s.coord(:,1) = X(isD);
%             s.coord(:,2) = Y(isD);
%             s.connec = delaunay(s.coord); 
%             m = Mesh(s);            
%             quad = Quadrature.set(m.type);
%             quad.computeQuadrature('CONSTANT');
%             dV = m.computeDvolume(quad);
%             isElemIns = dV < 1.0001*h^2/2;
%             s.connec = m.connec(isElemIns,:);
%             s.coord  = m.coord;
%             m2 = Mesh(s);
%             m = m2;
%             m.plot
             h = 0.03;
             model = createpde;
             vertCoord = [-1 1; 1 1; 1 0; -1 0];
             x = [-1,1,1,-1];
             y = [1,1,0,0];
             R1 = [3,4,vertCoord(:)']';

             centX = 0;
             centY = 0;
             rad = 0.1;
             C1 = [1,centX,centY,rad]';
             C1 = [C1;zeros(length(R1) - length(C1),1)];
             gm = [R1,C1];
             sf = 'R1-C1';
             ns = char('R1','C1');
             ns = ns';
             g = decsg(gm,sf,ns);

             geometryFromEdges(model,g);
             me = generateMesh(model,'Hmin',h,'Hmax',2*h,'GeometricOrder','linear');

             pdegplot(model,'EdgeLabels','on')
             axis equal
             xlim([-1.1,1.1])
             pdeplot(me)

             
             s.coord  = me.Nodes';
             s.connec = me.Elements';
             m = Mesh(s);
             m.plot()
        end

        function createBoundaryMesh(obj)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            b = boundary(x,y,1);
            obj.boundaryMesh = b;
        end

        function storeOrientationAngle(obj)
%             d = obj.experimentData;
%             alpha0  = d.dataRes.AlphaGauss;
%           %  obj.plotAlpha0(alpha0);
%             alpha(:,1) = obj.interpolateOrientationAngle(alpha0(:,1));
%             alpha(:,2) = obj.interpolateOrientationAngle(alpha0(:,2));
%            
 
%             
          %  theta = atan2(x2,x1);            
           
     %       thetaV = pi/6;
     %       theta = thetaV*ones(size(x2));

%      % Lineal + noise
%              x1 = obj.mesh.coord(:,1);    
%              x2 = obj.mesh.coord(:,2);
%              theta = 2*rand(1)*x1 + 2*rand(1)*x2 + rand(1) + 0.5*(2*rand(size(x1))-1)/2;
     
     % Pi discontinoity + noise
%              x1 = obj.mesh.coord(:,1);
%              x2 = obj.mesh.coord(:,2);
%              theta = zeros(size(x1)) + (pi/2);
%              theta(x1<0.5) = pi + (pi/2);
%              theta = theta + 0.5*(2*rand(size(x1))-1)/2;
%            
       %%%
     % Radial + noise
             x1 = obj.mesh.coord(:,1);
             x2 = obj.mesh.coord(:,2);
             theta = atan2(x2,x1) + 1*(2*rand(size(x1))-1)/2;
%%%

            alpha(:,1) = cos(theta);
            alpha(:,2) = sin(theta);            
            
            
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
            q = quiver(x,y,tx,ty);
            q.ShowArrowHead = 'off';
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
            q = quiver(x,y,tx,ty);
            q.ShowArrowHead = 'off';
            figure(2)
            subplot(1,2,iFigure)
            s.mesh = obj.mesh;
            s.field = atan2(t(:,1),t(:,2));
            p = NodalFieldPlotter(s);
            p.plot()
        end

        function project(obj)
            aBar = obj.orientationVector;            

            figure(23)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            aXBar = aBar(:,1);
            aYBar = aBar(:,2);
            q = quiver(x,y,aXBar,aYBar);
            q.ShowArrowHead = 'off';
            alpha = atan2(aYBar,aXBar);
            beta  = 2*alpha;
            bBar(:,1) = cos(beta);
            bBar(:,2) = sin(beta);


            figure(24)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = bBar(:,1);
            ty = bBar(:,2);
            quiver(x,y,tx,ty)


            b      = bBar;
            lambda = obj.computeInitialLambda();
            isErrorLarge = true;
            i = 1;
            while isErrorLarge
                cost(i)      = obj.computeCost(b,bBar);
              %  optPrimal(i) = obj.computePrimalOptimaility(lambda,u,alpha0);
              %  optDual(i)   = obj.computeDualHarmonicOptimality(u);
              %  error = norm([optPrimal,optDual]);
                
%                 theta = atan2(u(:,2),u(:,1));
%                 u2(:,1) = cos(2*theta);
%                 u2(:,2) = sin(2*theta);
                [uNew,lambda] = obj.solveProblem(bBar,b);
                beta = atan2(uNew(:,2),uNew(:,1));                               
                uNew = obj.projectInUnitBall(uNew);
                
             %   uNew2 = u2;
                
                err(i) = norm(uNew(:)-b(:))/norm(b(:));
                
                b = uNew;
                figure(200)
                clf
                plot(cost,'-+')                

                figure(201)
                plot(err,'-+')    

                isErrorLarge = err(i) > 1e-13;
                i = i +1;
                obj.plotOrientation(b,2) 
                obj.orientationAngle = beta;

            end

            
            bOpt = b;
            obj.plotOrientationBarAndOptimal(bBar,bOpt);            

            alpha = beta/2;
            aOpt(:,1) = cos(alpha);
            aOpt(:,2) = sin(alpha);
            obj.plotOrientationBarAndOptimal(aBar,aOpt);
           
        end

        function plotOrientationBarAndOptimal(obj,bBar,bOpt)
            figure()
            subplot(1,2,1)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = bBar(:,1);
            ty = bBar(:,2);
            q = quiver(x,y,tx,ty);
         %   q.ShowArrowHead = 'off';
            subplot(1,2,2)
            tx = bOpt(:,1);
            ty = bOpt(:,2);            
            q = quiver(x,y,tx,ty);
         %   q.ShowArrowHead = 'off';


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

       function dehomogenize(obj)
            obj.createBackgroundMesh();
            obj.createSuperEllipseParams();
            s.backgroundMesh     = obj.backgroundMesh;
            s.nCells             = 46;
            s.theta              = obj.orientationAngle;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
       end        

        function s = createLevelSetCellParams(obj)
            s.type   = 'smoothRectangle';
            s.widthH = obj.superEllipse.m1;
            s.widthV = obj.superEllipse.m2;
            s.pnorm  = obj.superEllipse.q;
            s.ndim   = 2;
        end       

        function createBackgroundMesh(obj)
            FV.vertices = [obj.mesh.coord,zeros(size(obj.mesh.coord,1),1)];
            FV.faces    = obj.mesh.connec;
            FV2 = refinepatch(FV);
            FV2 = refinepatch(FV2);
            s.coord = FV2.vertices(:,1:2);
            s.connec = FV2.faces;
            m = Mesh(s);
            obj.backgroundMesh = m;
        end        

        function createSuperEllipseParams(obj)
            sE.m1 = obj.interpolateFunction(obj.experimentData.dataRes.DesignVar1);
            sE.m2 = obj.interpolateFunction(obj.experimentData.dataRes.DesignVar2);
            sE.q  = obj.interpolateFunction(obj.experimentData.dataRes.SuperEllipseExponent);
            obj.superEllipse = sE;
        end
        
        function vq = interpolateFunction(obj,v)
            X = obj.mesh.coord(:,1);
            Y = obj.mesh.coord(:,2);
            F = scatteredInterpolant(X,Y,v);
            xB = obj.backgroundMesh.coord(:,1);
            yB = obj.backgroundMesh.coord(:,2);
            vq = F(xB,yB);
        end        

    end

end