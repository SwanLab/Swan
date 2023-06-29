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
          %  obj.trySomething();
            obj.createHarmonicProjection();
           % obj.createUnitBallProjector();
           % obj.storeOrientationAngle();
          %  obj.dehomogenize(); 
           % obj.computeSingularities();
            obj.project()
            obj.computeSingularities();
            obj.dehomogenize();                        
        end

    end

    methods (Access = private)

        function init(obj)
            close all
       %     obj.filePath = 'Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
       %     obj.fileName = 'LshapeCoarseSuperEllipseDesignVariable';
       %     obj.iteration = 665;
% 
            obj.filePath = 'Topology Optimization/Applications/Dehomogenizing/ExampleCompliance/';  
             obj.fileName = 'ExperimentingPlotSuperEllipse';
             obj.iteration = 64;
                        
        end

        function loadDataExperiment(obj)
            s.fileName = [obj.fileName,num2str(obj.iteration)];
            s.folderPath = fullfile(obj.filePath );
            w = WrapperMshResFiles(s);
            w.compute();
            obj.experimentData = w;
        end

        function trySomething(obj)
            
            %o = obj.tryOrientationCase1(sigma0);

         


            oS2 = obj.tryOrientationCase2(sigma1);
  



            %a1d(:,1,:) = obj.experimentData.dataRes.AlphaGauss';
            %tFd = obj.createOrientation(a1d);
          

        end

        function tF = createVectorFromSigma(obj)
            s.mesh    = obj.experimentData.mesh;
            s.fValues = obj.experimentData.dataRes.StressPrimal;
            sigma0 = P0Function(s);
            sigma1 = sigma0.project('P1');
            

            s.type = '2D';
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            s2(1,:,:) = (sigma1.fValues');
            pcomp.compute(s2);

            a = pcomp.direction;
            a1 = a(:,1,:);
            a2 = a(:,2,:);

            s.fValues = squeeze(a1)';
            s.mesh    = obj.mesh;
            aF{1} = P1Function(s);

            s.fValues = squeeze(a2)';
            s.mesh    = obj.mesh;
            aF{2} = P1Function(s);
            tF = aF;
        end

        function o = tryOrientationCase2(obj,sigma1)

            tF = obj.createVectorFromSigma(sigma1)
         %   x1 = (a1(1,:));
         %   x2 = (a1(2,:));
         %   tV = atan2(x2,x1);
         %   s.mesh = obj.mesh;
         %   s.fValues = squeeze(tV');
         %   tF = P1Function(s);
         %   tF = tF.project('P0');
            

           % tF = obj.createOrientation(a1);
            
            s.theta = tF;%
            s.mesh  = obj.mesh;
            o = OrientationVectors(s);   
          %  o.computeDeformedCoordinates 
        end  


        function o = tryOrientationCase1(obj,sigma0)
            s.type = '2D';
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            s1 = permute(sigma0.fValues,[2 1 3]);
            pcomp.compute(s1);

            a = pcomp.direction;
            a1(:,1,:) = a(:,1,:);

            
            tF = obj.createOrientation(a1);
            
            s.theta = tF;
            s.mesh  = obj.mesh;
            o = OrientationVectors(s);   
          %  o.computeDeformedCoordinates 
        end

        function createMesh(obj)
            d = obj.experimentData;
            obj.mesh = d.mesh;          
         % obj.mesh = obj.createSquareMesh();
         %  obj.mesh  = obj.createRectangularWithCircularHoleMesh();
        end

      function aBar = createOrientationByHand(obj)
           % alpha = obj.createLinealWithNoiseAngle();
            alpha = obj.createRadialWithNoiseAngle();
           % alpha = obj.createPiDiscontinuityAndNoiseAngle();
            aBar(:,1) = cos(alpha);
            aBar(:,2) = sin(alpha); 
      end

        function tF = createOrientation(obj,a)
            x1 = (a(1,:,:));
            x2 = (a(2,:,:));
            tV = atan2(x2,x1);
            s.mesh = obj.mesh;
            s.fValues = squeeze(tV);
            tF = P0Function(s);
        end      
      

        function [a1,b1] = createOrientationVector(obj)
         %   aBar = obj.createOrientationByHand();
            a0 = obj.createOrientationFromData();

            s.theta = obj.createOrientation(a0.fValues);
            s.mesh  = obj.mesh;
            o = OrientationVectors(s);            

            a0X = squeeze(a0.fValues(1,1,:));
            a0Y = squeeze(a0.fValues(2,1,:));
            b0V = obj.createDobleOrientationVector(a0X,a0Y);
            s.fValues = b0V;
            s.mesh = obj.mesh;
            b0 = P0Function(s);
            b1a = b0.project('P1');

            a1   = a0.project('P1');  
            a1X = a1.fValues(:,1);
            a1Y = a1.fValues(:,2);
            b1V = obj.createDobleOrientationVector(a1X,a1Y);
            s.fValues = b1V;
            s.mesh = obj.mesh;
            b1b = P1Function(s);            

            
            b1 = b1b;
            
         %   b = obj.computeHalfOrientationFromOrientation(a);
          %  aBar = obj.interpolateOrientation(aBar);
          %  bBar = obj.interpolateOrientation(bBar);
        end        

        function m = createSquareMesh(obj)
            h = 0.05;
            [X,Y] = meshgrid(-1:h:1,0:h:1); 
           s.coord(:,1) = X(:);
           s.coord(:,2) = Y(:);
           s.connec = delaunay(s.coord); 

%             [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
%             s.coord  = V(:,1:2);
%               s.connec = F;

            m = Mesh(s);            
            m.plot;
        end

        function m = createRectangularWithCircularHoleMesh(obj)
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

        % function createOrientationVector(obj)
        %     t       = obj.theta.fValues;
        %     a1(:,1) = cos(t);
        %     a1(:,2) = sin(t);
        %     a2(:,1) = -sin(t);
        %     a2(:,2) = cos(t);
        %     a(:,:,1) = a1;
        %     a(:,:,2) = a2;
        %     nDim = obj.mesh.ndim;
        %     orientation = cell(nDim,1);
        %     for iDim = 1:nDim
        %         s.fValues = a(:,:,iDim);
        %         s.mesh   = obj.mesh;
        %         bf = P1Function(s);
        %         orientation{iDim} = bf;
        %     end 
        %     obj.orientationVectorFunc = orientation;
        % end

        function alpha = createLinealWithNoiseAngle(obj)
             coord0 = obj.mesh.computeBaricenter';
             x1 = coord0(:,1);
             x2 = coord0(:,2);
             alpha = 2*rand(1)*x1 + 2*rand(1)*x2 + rand(1) + 0.5*(2*rand(size(x1))-1)/2;
        end

        function alpha = createRadialWithNoiseAngle(obj)
             coord0 = obj.mesh.computeBaricenter';
             x1 = coord0(:,1);
             x2 = coord0(:,2);
             alpha = atan2(x2+0.5,x1) + 0*(2*rand(size(x1))-1)/2;
        end

        function alpha = createPiDiscontinuityAndNoiseAngle(obj)
             coord0 = obj.mesh.computeBaricenter';
             x1 = coord0(:,1);
             x2 = coord0(:,2);
             alpha = zeros(size(x1)) + (pi/2);
             alpha(x1<0) = pi + (pi/2);
             alpha = alpha;% + 0.5*(2*rand(size(x1))-1)/2;
        end

        function aF = createOrientationFromData(obj)
            d = obj.experimentData;
            aBar  = d.dataRes.AlphaGauss;
            s.mesh    = obj.mesh;
            s.fValues = aBar;
            aF = P0Function(s);
        end        

        function storeOrientationAngle(obj)
            [aBar,bBar] = obj.createOrientationVector();            
            aBar = obj.projectInUnitBall(aBar);
            alpha(:,1) = atan2(aBar(:,1),aBar(:,2));  
            obj.plotOrientation(aBar,1);
            obj.orientationAngle  = alpha;
            obj.orientationVector = aBar;            
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

        function plotOrientationVector(obj,b)
            a = obj.createHalfOrientationVector(b(:,1),b(:,2));
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = a(:,1);
            ty = a(:,2);
            q = quiver(x,y,tx,ty);
            q.ShowArrowHead = 'off';            
        end

        % function a1 = interpolateOrientation(obj,a0)
        %     a1(:,1) = obj.interpolateOrientationComponent(a0(:,1));
        %     a1(:,2) = obj.interpolateOrientationComponent(a0(:,2));
        % end
        % 
        % function vI = interpolateOrientationComponent(obj,v0)
        %     s.mesh    = obj.mesh;
        %     s.fValues = v0;
        %     p = PieceWiseConstantFunction(s);
        %     vI = p.projectToLinearNodalFunction(); 
        % end


        function plotOrientation(obj,t,iFigure)
            figure(1)
            subplot(1,2,iFigure)
            obj.plotOrientationVector(t);
            figure(2)
            subplot(1,2,iFigure)  
            obj.plotAngleField(t);
        end

        function plotAngleField(obj,t)
            s.mesh = obj.mesh;
            s.field = atan2(t(:,2),t(:,1));
            p = NodalFieldPlotter(s);
            p.plot()
            caxis([-pi pi])
            colormap default
            cmap_pos = colormap;
            cmap_neg = flipud(cmap_pos);
            cmap = [cmap_neg; cmap_pos];
            colormap(cmap)
            t= {'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'};
            c = colorbar;
            c.Ticks = -pi:pi/4:pi;
            c.TickLabels = t;
        end


        function b = createDobleOrientationVector(obj,aX,aY)
            alpha = atan2(aY,aX);
            beta  = 2*alpha;
            b(:,1) = cos(beta);
            b(:,2) = sin(beta);
        end

        function a = createHalfOrientationVector(obj,bX,bY)
            beta   = atan2(bY,bX);
            alpha  = beta/2;
            a(:,1) = cos(alpha);
            a(:,2) = sin(alpha);
        end    

        function rho = computeRho(obj)            
            %rho = obj.experimentData.dataRes.DensityGauss;           
            %rho = obj.interpolateOrientationComponent(rho);
            rho = zeros(size(obj.mesh.coord(:,1)));
            rho(:) = 0.5;
        end

        function project(obj)
            aBar = obj.createVectorFromSigma();
            a1 = aBar{1};
            a1X = a1.fValues(:,1);
            a1Y = a1.fValues(:,2);
            bBar = obj.createDobleOrientationVector(a1X,a1Y);
           
            s.fValues = bBar;
            s.mesh    = obj.mesh;
            bBarF = P1Function(s);
          %  bBarF.plotArrowVector();


            %[aBar,bBar] = obj.createOrientationVector(); 

            figure(23)
            obj.plotOrientationVector(bBar);

            
            rho = obj.computeRho();

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
                [uNew,lambda] = obj.solveProblem(rho,bBar,b);
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

                isErrorLarge = err(i) > 1e-2;%1e-13;
                i = i +1;
                obj.plotOrientation(b,2) 
                alpha = beta/2;
                obj.orientationAngle = alpha;

            end

            
            bOpt = b;
            obj.plotOrientationBarAndOptimal(bBar,bOpt);            

            alpha = beta/2;
            aOpt(:,1) = cos(alpha);
            aOpt(:,2) = sin(alpha);
            obj.plotOrientationBarAndOptimal(aBar,aOpt);
           
            figure()
            subplot(1,2,1)
            obj.plotAngleField(bBar)
            title('$\bar{\beta}$','Interpreter','latex','FontSize',15)
            subplot(1,2,2)
            obj.plotAngleField(b)
            title('${\beta}^*$','Interpreter','latex','FontSize',15)            
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

        function [v,lambda] = solveProblem(obj,rho,alpha0,vH)
           h  = obj.harmonicProjector;
           [v,lambda] = h.solveProblem(rho,alpha0,vH);
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

        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation(:,1) = cos(obj.orientationAngle);
            s.orientation(:,2) = sin(obj.orientationAngle);
            sF = SingularitiesFinder(s);
            isS = sF.computeSingularElements();
            sF.plot();
        end        

       function dehomogenize(obj)
            obj.createBackgroundMesh();
         %   obj.createSuperEllipseParams();
            s.backgroundMesh     = obj.backgroundMesh;
            s.nCells             = 46;

%              x1 = obj.backgroundMesh.coord(:,1);
%              x2 = obj.backgroundMesh.coord(:,2);
%              alpha = zeros(size(x1)) + (pi/2);
%              alpha(x1<0) = pi + (pi/2);
%              alpha = alpha;
% 
%             
             alpha  = obj.createRadialWithNoiseAngle();                
             alpha = obj.interpolateOrientationComponent(alpha);
             x1 = obj.mesh.coord(:,1);
             x2 = obj.mesh.coord(:,2);   
             alpha(x1<0) = alpha(x1<0)  + (pi);

             s.theta              = alpha;
%            
%                 v = 2*obj.orientationAngle;
%                 X = obj.mesh.coord(:,1);
%                 Y = obj.mesh.coord(:,2);
%                 F = scatteredInterpolant(X,Y,v);
%                 xB = obj.backgroundMesh.coord(:,1);
%                 yB = obj.backgroundMesh.coord(:,2);
%                 vq = F(xB,yB);
%                 obj.mesh = obj.backgroundMesh;
%                 s.theta = vq/2;

            

   %         s.theta              = obj.orientationAngle;
           
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
       end        

        function s = createLevelSetCellParams(obj)
           % s.type   = 'smoothRectangle';
            s.type   = 'rectangleInclusion';
            s.widthH = 0.7*ones(size(obj.backgroundMesh.coord,1),1);
            s.widthV = 0.7*ones(size(obj.backgroundMesh.coord,1),1);
            %s.pnorm  = obj.superEllipse.q;
            s.ndim   = 2;
        end       

        function createBackgroundMesh(obj)
            FV.vertices = [obj.mesh.coord,zeros(size(obj.mesh.coord,1),1)];
            FV.faces    = obj.mesh.connec;
            FV2 = FV;
            FV2 = refinepatch(FV2);
            FV2 = refinepatch(FV2);
            FV2 = refinepatch(FV2);
          %  FV2 = refinepatch(FV2);
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