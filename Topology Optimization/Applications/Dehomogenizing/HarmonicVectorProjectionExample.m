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
        orientationVector
        doubleOrientationVector        
        harmonicProjector
        unitBallProjector
    end
    
    methods (Access = public)

        function obj = HarmonicVectorProjectionExample()
            obj.init();
            obj.loadDataExperiment();
            obj.createMesh();
            obj.createBoundaryMesh();
            obj.createOrientationVector();
          %  obj.createDoubleOrientationVector();
            obj.createHarmonicProjection();
            obj.createUnitBallProjector();
            obj.project()
            obj.computeSingularities();
            obj.dehomogenize();                        
        end

    end

    methods (Access = private)

        function init(obj)
            close all
  %          obj.filePath = 'Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
  %          obj.fileName = 'LshapeCoarseSuperEllipseDesignVariable';
  %          obj.iteration = 665;
% 
             obj.filePath = 'Topology Optimization/Applications/Dehomogenizing/ExampleCompliance/';  
             obj.fileName = 'ExperimentingPlotSuperEllipse';
             obj.iteration = 64;
        end

        function loadDataExperiment(obj)
           % s.fileName = [obj.fileName,num2str(obj.iteration)];
           %  s.folderPath = fullfile(obj.filePath );
           %  w = WrapperMshResFiles(s);
           %  w.compute();
            d = load('DataExample.mat');
            w = d.w;
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

        function createOrientationVector(obj)
         %   a1 = obj.createOrientationVectorFromProjectedSigma();            
            a0 = obj.createOrientationFromSigma();
            obj.orientationVector = a0;
        end

        function a = createOrientationFromSigma(obj)
            sigma0  = obj.getSigma0FromData();
            sigma0V(1,:,:) = sigma0.fValues';
            [a1,a2] = obj.obtainPrincipalDirections(sigma0V);   
            a{1}    = obj.createFunctionP0(a1);
            a{2}    = obj.createFunctionP0(a2);
        end              

       function a  = createOrientationVectorFromProjectedSigma(obj)
            sigma0  = obj.getSigma0FromData();            
            sigma1 = sigma0.project('P1');
            sigma1V(1,:,:) = sigma1.fValues';
            [a1,a2] = obj.obtainPrincipalDirections(sigma1V);
            a{1}    = obj.createFunctionP1(a1);
            a{2}    = obj.createFunctionP1(a2);            
        end        

        function sigma0 = getSigma0FromData(obj)
            s.mesh    = obj.experimentData.mesh;
            s.fValues = obj.experimentData.dataRes.StressPrimal;
            sigma0 = P0Function(s);
        end

        function [a1,a2] = obtainPrincipalDirections(obj,sigma)
            s.type = '2D';
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            pcomp.compute(sigma);
            a = pcomp.direction;
            a1 = a(:,1,:);
            a2 = a(:,2,:);
        end

        function aF = createFunctionP0(obj,a)
            s.mesh = obj.mesh;
            s.fValues = squeeze(a)';
            aF = P0Function(s);
        end                     

        function aF = createFunctionP1(obj,a)
            s.fValues = squeeze(a);
            s.mesh    = obj.mesh;
            aF = P1Function(s);
        end

        function createHarmonicProjection(obj)
            s.mesh         = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;            
            h = LinearizedHarmonicProjector2(s);
            obj.harmonicProjector = h;
        end     

        function aF = createOrientationFromData(obj)
            d = obj.experimentData;
            aBar  = d.dataRes.AlphaGauss;
            s.mesh    = obj.mesh;
            s.fValues = aBar;
            aF = P0Function(s);
        end        

        function plotOrientationVector(obj,b)
            a = obj.createHalfOrientationVectorP0(b(:,1),b(:,2));
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = a(:,1);
            ty = a(:,2);
            q = quiver(x,y,tx,ty);
            q.ShowArrowHead = 'off';            
        end


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


        function b = createDobleOrientationVectorP1(obj,a)
            aX = squeeze(a.fValues(:,1));
            aY = squeeze(a.fValues(:,2));            
            alpha = atan2(aY,aX);
            beta  = 2*alpha;
            bV(:,1) = cos(beta);
            bV(:,2) = sin(beta);       
            s.fValues = bV;
            s.mesh    = obj.mesh;
            b = P1Function(s);
        end

        function b0 = createDoubleOrientationVectorP0(obj,a0)
            aX = squeeze(a0.fValues(1,1,:));
            aY = squeeze(a0.fValues(2,1,:));
            alpha = atan2(aY,aX);
            beta  = 2*alpha;
            b(:,1) = cos(beta);
            b(:,2) = sin(beta);                        
            s.fValues = b;
            s.mesh    = obj.mesh;
            b0 = P0Function(s);
        end

        function a = createHalfOrientationVectorP0(obj,b)
            bX = squeeze(b.fValues(1,1,:));
            bY = squeeze(b.fValues(2,1,:));
            beta   = atan2(bY,bX);
            alpha  = beta/2;
            aV(:,1) = cos(alpha);
            aV(:,2) = sin(alpha);
            s.fValues = aV;
            s.mesh    = obj.mesh;
            a = P0Function(s);
        end    

        function rho = computeRho(obj)            
            %rho = obj.experimentData.dataRes.DensityGauss;           
            %rho = obj.interpolateOrientationComponent(rho);
            rhoV = zeros(size(obj.mesh.coord(:,1)));
            rhoV(:) = 0.5;
            s.fValues = rhoV;
            s.mesh    = obj.mesh;
            rho = P1Function(s);
        end

        function b1 = projectViaFilterIteration(obj,rho,b1)
            bInitial = b1;

            h    = obj.harmonicProjector;
            [hRes,lRes] = h.evaluateAllResiduals(rho,bInitial,b1);
            hRes.plot();
            colorbar
            hMax = max(hRes.fValues(:));
            hMin = min(hRes.fValues(:));

            lRes.plot();
            colorbar


            a1 = obj.createHalfOrientationVectorP1(b1);
            for i = 1:5000
              %  figure(1)
              %  a1.plot();
              %  figure(2)
              %  a1.plotArrowVector();
                b1 = b1.project('H1P1');

                b1 = obj.projectInUnitBall(b1.fValues);
                b1 = obj.createFunctionP1(b1);

                a1 = obj.createHalfOrientationVectorP1(b1);
                a1 = obj.projectInUnitBall(a1.fValues);
                a1 = obj.createFunctionP1(a1);               

                s.mesh        = obj.mesh;
                s.orientation = a1;
                sC = SingularitiesComputer(s);
                sC.compute();

                

                            

                if ~sC.isSingularityInBoundary() && ~sC.isNodeInTwoSingularElements() && sC.nSing < 1
                    break
                    sC.plot
                    %%% NOT Touching!!! 
                   % a = 0;

                    %  a{1} = a0;
                    % 
                    %  s.fValues(:,2) = a0.fValues(1,1,:);
                    %  s.fValues(:,1) = -a0.fValues(2,1,:);
                    %  s.mesh           = obj.mesh;
                    %  a{2}             = P0Function(s);
                    % obj.dehomogenize(a);
                end
            end
        end

        function bBar = projectViaPicard(obj,rho,bBar,bInitial)  
            h  = obj.harmonicProjector;
            [bBar,lambda] = h.solveProblem(rho,bBar,bInitial);
        end

        function a1 = createHalfOrientationVectorP1(obj,b1)
            bX = b1.fValues(:,1);
            bY = b1.fValues(:,2);
            beta   = atan2(bY,bX);
            alpha  = beta/2;
            aV(:,1) = cos(alpha);
            aV(:,2) = sin(alpha);
            s.fValues = aV;
            s.mesh    = obj.mesh;
            a1 = P1Function(s);
        end

        function project(obj)
            a = obj.orientationVector();
            a0 = a{1};
            a1 = a0.project('P1');
            a1 = obj.projectInUnitBall(a1.fValues);
            a1 = obj.createFunctionP1(a1);            


            b1 = obj.createDobleOrientationVectorP1(a1);
            b1 = obj.projectInUnitBall(b1.fValues);
            b1 = obj.createFunctionP1(b1);            
        
            figure(3)
            s.mesh        = obj.mesh;
            s.orientation = a1;
            sC = SingularitiesComputer(s);
            sC.compute();
            sC.plot();         

            rho = obj.computeRho();                        
                  
            h = obj.harmonicProjector;
            [hRes,lRes] = h.evaluateAllResiduals(b1,b1);
            hRes.plot();
            lRes.plot();                 
        

            %a01 = obj.createHalfOrientationVector(b0);
      
            %%%Via filtering 
            bNew = obj.projectViaFilterIteration(rho,b1);
       %     a1   = obj.createHalfOrientationVectorP1(bNew);
            
%              h    = obj.harmonicProjector;
%              [hRes,lRes] = h.evaluateAllResiduals(rho,b1,bNew);
%              hRes.plot();
%              lRes.plot();

            
            
            %%%Via picard
            %b1   = bNew;
      %      bNew = obj.projectViaPicard(rho,b1,b1);
            a1   = obj.createHalfOrientationVectorP1(bNew);
            
            h    = obj.harmonicProjector;
            hRes = h.evaluateAllResiduals(rho,b1,bNew);
            hRes.plot();



            
            a1 = obj.projectInUnitBall(a1.fValues);
            a1 = obj.createFunctionP1(a1);            

            figure(4)
            s.mesh        = obj.mesh;
            s.orientation = a1;
            sC = SingularitiesComputer(s);
            sC.compute();
            sC.plot();


            a{1} = a1;

            s.fValues(:,2) = a1.fValues(:,1);
            s.fValues(:,1) = -a1.fValues(:,2);
            s.mesh           = obj.mesh;
            a{2}             = P1Function(s);
            obj.dehomogenize(a);


            b1 = obj.doubleOrientationVector{1};
            b1.plotArrowVector();



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

       function dehomogenize(obj,a1)
%             obj.createBackgroundMesh();
%          %   obj.createSuperEllipseParams();
%             s.backgroundMesh     = obj.backgroundMesh;
%             s.nCells             = 46;
% 
% %              x1 = obj.backgroundMesh.coord(:,1);
% %              x2 = obj.backgroundMesh.coord(:,2);
% %              alpha = zeros(size(x1)) + (pi/2);
% %              alpha(x1<0) = pi + (pi/2);
% %              alpha = alpha;
% % 
% %             
%              alpha  = obj.createRadialWithNoiseAngle();                
%              alpha = obj.interpolateOrientationComponent(alpha);
%              x1 = obj.mesh.coord(:,1);
%              x2 = obj.mesh.coord(:,2);   
%              alpha(x1<0) = alpha(x1<0)  + (pi);
% 
%              s.theta              = alpha;
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
            s.nCells  = 55;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.theta              = a1;
            s.mesh               = obj.mesh;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
       end        

        function s = createLevelSetCellParams(obj)

        

            s.type   = 'smoothRectangle';
           % s.type   = 'rectangleInclusion';
            s.widthH = obj.experimentData.dataRes.DesignVar1;%0.85*ones(size(obj.mesh.coord,1),1);
            s.widthV = obj.experimentData.dataRes.DesignVar2;%0.85*ones(size(obj.mesh.coord,1),1);
            s.pnorm  = obj.experimentData.dataRes.SuperEllipseExponent;
            s.ndim   = 2;
        end       








%%%%%%%%%%%%%%%%%%%%%%%%


    function aBar = createOrientationByHand(obj)
           % alpha = obj.createLinealWithNoiseAngle();
            alpha = obj.createRadialWithNoiseAngle();
           % alpha = obj.createPiDiscontinuityAndNoiseAngle();
            aBar(:,1) = cos(alpha);
            aBar(:,2) = sin(alpha); 
      end

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


   
    end

end