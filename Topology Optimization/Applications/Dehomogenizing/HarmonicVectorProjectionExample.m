classdef HarmonicVectorProjectionExample < handle

    properties (Access = public)

    end

    properties (Access = private)
        harmonicVector
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
             obj.harmonize();
            %obj.harmonizeWithPenalizedUnitNorm();
            %obj.harmonizeWithPenalizedHarmonizity();
          
            obj.dehomogenize();                        
        end

    end

    methods (Access = private)

        function init(obj)
            close all
           % obj.filePath = 'Topology Optimization/Applications/Dehomogenizing/ExampleLShape/';
           %  obj.fileName = 'LshapeCoarseSuperEllipseDesignVariable';
           % obj.iteration = 665;
 
           obj.filePath = 'Topology Optimization/Applications/Dehomogenizing/ExampleCompliance/';  
           obj.fileName = 'ExperimentingPlotSuperEllipse';
           obj.iteration = 64;
        end

        function loadDataExperiment(obj)
     %      s.fileName = [obj.fileName,num2str(obj.iteration)];
     %      s.folderPath = fullfile(obj.filePath );
     %      w = WrapperMshResFiles(s);
     %      w.compute();
           d = load('DataExample.mat');
           w = d.w;
             obj.experimentData = w;
        end

        function createMesh(obj)
            d = obj.experimentData;
            s.coord  = d.mesh.coord;          
            s.connec = d.mesh.connec;         
            obj.mesh = Mesh.create(s);
        end    

        function createBoundaryMesh(obj)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            b = boundary(x,y,1);
            obj.boundaryMesh = b;
        end        

        function createOrientationVector(obj)
            sigma0  = obj.getSigma0FromData();
            sigma0V(1,:,:) = sigma0.fValues';
            [a1,a2] = obj.obtainPrincipalDirections(sigma0V);   
            a{1} = obj.createP1DiscontinousOrientation(a1);
            a{1} = obj.projectInUnitBall(a{1});
            a{2} = obj.createP1DiscontinousOrientation(a2);
            a{2} = obj.projectInUnitBall(a{2}); 
            obj.orientationVector = a;
        end
             

        function sigma0 = getSigma0FromData(obj)
            s.mesh    = obj.mesh;
            s.fValues = obj.experimentData.dataRes.StressPrimal;
            s.order   = 'P0';
            sigma0 = LagrangianFunction(s);
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

        function aF = createP1DiscontinousOrientation(obj,a)
            s.mesh = obj.mesh;
            s.fValues = squeeze(a)';
            s.order   = 'P0';
            aF = LagrangianFunction(s);
            aF = project(aF,'P1D');
        end                     

        function aF = createFunctionP1(obj,a)
            s.fValues = squeeze(a);
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            aF = LagrangianFunction(s);
        end

        function h = createHarmonicProjectionWithPenalizedUnitNorm(obj,epsilon)
            s.mesh         = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;  
            s.epsilon      = epsilon;
            h = LinearizedHarmonicProjector2(s);
        end             

        function h = createHarmonicProjection(obj)
            s.mesh         = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;  
            h = LinearizedHarmonicProjector3(s);
        end                   

        function h = createHarmonicProjectionWithPenalizedHarmonizity(obj,epsilon)
            s.mesh         = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;  
            s.epsilon      = epsilon;
            h = LinearizedHarmonicProjector4(s);
        end            




        function aF = createOrientationFromData(obj)
            d = obj.experimentData;
            aBar  = d.dataRes.AlphaGauss;
            s.mesh    = obj.mesh;
            s.order   = 'P0';
            aF = LagrangianFunction(s);
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
            s.order   = 'P1D';
            b = LagrangianFunction(s);
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
            s.order   = 'P0';
            a = LagrangianFunction(s);
        end    

        function rho = computeRho(obj)            
            %rho = obj.experimentData.dataRes.DensityGauss;           
            %rho = obj.interpolateOrientationComponent(rho);
            rhoV = zeros(size(obj.mesh.coord(:,1)));
            rhoV(:) = 0.5;
            s.fValues = rhoV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            rho = LagrangianFunction(s);
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
            s.order   = 'P0';
            a1 = LagrangianFunction(s);
        end

        function resHNorm = plotAll(obj,h,bBar,b)
            [resL,resH,resB,resG] = h.evaluateAllResiduals(bBar,b);

            a1   = obj.createHalfOrientationVectorP1(b);
            a1 = obj.projectInUnitBall(a1.fValues);
            a1 = obj.createFunctionP1(a1);
            figure()
            s.mesh        = obj.mesh;
            s.orientation = a1;
            sC = SingularitiesComputer(s);
            sC.compute();
            sC.plot();
            title('Singularities')
            resL.plot
            title('L2Distance')
            resH.plot
            title('Harmonicity')
            resB.plot
            title('UnitBall')
            resG.plot
            title('Gradient')
            resHNorm = resH.computeL2norm();
        end

        function harmonizeWithPenalizedHarmonizity(obj)
            bInit = obj.obtainInitialOrientationVector();        
            bBar  = bInit;

            epsilons = linspace(0,1,3);
           
            for i = 1:length(epsilons)
                epsilon = epsilons(i);
                
                for k = 1:1
                h = obj.createHarmonicProjectionWithPenalizedHarmonizity(epsilon);

                obj.plotAll(h,bBar,bInit);
                bNew = h.solveProblem(bBar,bInit);
                obj.plotAll(h,bBar,bNew);


                a1 = obj.createHalfOrientationVectorP1(bNew);
                a1 = obj.projectInUnitBall(a1.fValues);
                a1 = obj.createFunctionP1(a1);
                bNewP = obj.createDobleOrientationVectorP1(a1);



                hnorm(i) = obj.plotAll(h,bBar,bNewP);
                for j = 1:i
                disp(['epsilon ',num2str(epsilons(j)),' hNorm ',num2str(hnorm(j))])
                end
                close all

                bInit = bNew;
                end                

            end





            a{1} = a1;

            s.fValues(:,2) = a1.fValues(:,1);
            s.fValues(:,1) = -a1.fValues(:,2);
            s.mesh         = obj.mesh;
            s.order        = 'P1';
            a{2}           = LagrangianFunction(s);

            obj.harmonicVector = a;
            

        end

        function harmonizeWithPenalizedUnitNorm(obj)

            bInit = obj.obtainInitialOrientationVector();        
            bBar  = bInit;

            epsilons = linspace(0,1,3);
           
            for i = 1:length(epsilons)
                epsilon = epsilons(i);
                
                for k = 1:1
                h = obj.createHarmonicProjectionWithPenalizedUnitNorm(epsilon);

                obj.plotAll(h,bBar,bInit);
                bNew = h.solveProblem(bBar,bInit);
                obj.plotAll(h,bBar,bNew);


                a1 = obj.createHalfOrientationVectorP1(bNew);
                a1 = obj.projectInUnitBall(a1.fValues);
                a1 = obj.createFunctionP1(a1);
                bNewP = obj.createDobleOrientationVectorP1(a1);



                hnorm(i) = obj.plotAll(h,bBar,bNewP);
                for j = 1:i
                disp(['epsilon ',num2str(epsilons(j)),' hNorm ',num2str(hnorm(j))])
                end
                close all

                bInit = bNew;
                end                

            end



            a{1} = a1;

            s.fValues(:,2) = a1.fValues(:,1);
            s.fValues(:,1) = -a1.fValues(:,2);
            s.mesh         = obj.mesh;
            s.order        = 'P1';
            a{2}           = LagrangianFunction(s);

            obj.harmonicVector = a;
            

               
        end

        function harmonize(obj)

            a  = obj.orientationVector();
            a1 = a{1};
            bInit = obj.createDobleOrientationVectorP1(a1);             
            bBar  = bInit;

            h = obj.createHarmonicProjection();

               
            for k = 1:1

                obj.plotAll(h,bBar,bInit);
                bNew = h.solveProblem(bBar,bInit);
                obj.plotAll(h,bBar,bNew);


                a1 = obj.createHalfOrientationVectorP1(bNew);
                a1 = obj.projectInUnitBall(a1);
                bNewP = obj.createDobleOrientationVectorP1(a1);



                hnorm(k) = obj.plotAll(h,bBar,bNewP);
                for j = 1:k
                    disp(['kIter ',num2str(k),' hNorm ',num2str(hnorm(j))])
                end
                close all

                bInit = bNew;
            end




            a{1} = a1;

            s.fValues(:,2) = a1.fValues(:,1);
            s.fValues(:,1) = -a1.fValues(:,2);
            s.mesh         = obj.mesh;
            s.order        = 'P1';
            a{2}           = LagrangianFunction(s);

            obj.harmonicVector = a;




        end


        function [v,lambda] = solveProblem(obj,rho,alpha0,vH)
           h  = obj.harmonicProjector;
           [v,lambda] = h.solveProblem(rho,alpha0,vH);
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
            
          %  s = load('HarmonicVector');
            s = load('HarmonicVectorCantilever');
            obj.harmonicVector = s.a;
            obj.mesh = s.m;
            %obj.experimentData = s.e;

            a = obj.harmonicVector;
                
            s.nCells  = 75;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.theta              = a;
            s.mesh               = obj.mesh;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();


           s.nCells             = 75;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            s.orientationA       = a; 
            d = Dehomogenizer(s);
            ls = d.compute();
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


    methods (Access = private, Static)

        function vNF = projectInUnitBall(vF)
            v    = vF.fValues;  
            norm = vecnorm(v,2,2);
            v    = v./norm;
            vNF = LagrangianFunction.create(vF.mesh,vF.ndimf,vF.order);
            vNF.setFValues(v);
        end        

    end

end