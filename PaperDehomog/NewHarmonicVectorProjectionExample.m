classdef NewHarmonicVectorProjectionExample < handle

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

        function obj = NewHarmonicVectorProjectionExample()
            obj.init();
            obj.loadDataExperiment();
            obj.createMesh();
            obj.createBoundaryMesh();
            obj.createOrientationVector();
            obj.harmonize();
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
            s1 = project(sigma0,'P1');
            %sigma0V(1,:,:) = sigma0.fValues';
            a = obj.obtainPrincipalDirections(s1);
            obj.orientationVector = a;
        end

        function sigma0 = getSigma0FromData(obj)
            s.mesh    = obj.mesh;
            s.fValues = obj.experimentData.dataRes.StressPrimal;
            s.order   = 'P0';
            sigma0 = LagrangianFunction(s);
        end

        function a = obtainPrincipalDirections(obj,sigma)
            s.type = '2D';
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            [a,~] = pcomp.compute(sigma);
            a{1} = obj.projectInUnitBall(a{1});
            a{2} = obj.projectInUnitBall(a{2});
        end

        function h = createHarmonicProjection(obj)
            s.mesh         = obj.mesh;
            s.boundaryMesh = obj.boundaryMesh;
            h = LinearizedHarmonicProjector3(s);
        end

        function harmonize(obj)
            a     = obj.orientationVector;
            a1    = a{1};
            bInit = obj.createDobleOrientationVectorP1(a1);
            bBar  = bInit;
            h     = obj.createHarmonicProjection();

            obj.plotAll(h,bBar,bInit);
            bNew = h.solveProblem(bBar,bInit);
            bNew = obj.projectInUnitBall(bNew);
            obj.plotAll(h,bBar,bNew);
            a1   = obj.createHalfOrientationVectorP1(bNew);
            
            a1OrtV(:,1) = -a1.fValues(:,2);
            a1OrtV(:,2) = a1.fValues(:,1);
            a1Ort = LagrangianFunction.create(obj.mesh,2,'P1');
            a1Ort.setFValues(a1OrtV);
            obj.harmonicVector{1} = a1;
            obj.harmonicVector{2} = a1Ort;
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
            s.order   = 'P1';
            b = LagrangianFunction(s);
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
            s.order   = 'P1';
            a1 = LagrangianFunction(s);
        end

        function resHNorm = plotAll(obj,h,bBar,b)
            close all
            [resL,resH,resB,resG] = h.evaluateAllResiduals(bBar,b);

            a1   = obj.createHalfOrientationVectorP1(b);
            s.mesh        = obj.mesh;
            s.orientation = a1;
            sC = SingularitiesComputer(s);
            sCf = sC.compute();
            plot(sCf);
            title('Singularities')
            plot(resL)
            title('L2Distance')
            plot(resH)
            title('Harmonicity')
            plot(resB)
            title('UnitBall')
            plot(resG)
            title('Gradient')
            resHNorm = Norm(resH,'L2');

            figHandles = findall(groot, 'Type', 'figure');
            numFigures = length(figHandles);
            figure; tiledlayout(ceil(sqrt(numFigures)), ceil(sqrt(numFigures)));

            for i = 1:numFigures
                nexttile;
                axesHandles = findall(figHandles(i), 'Type', 'axes');
                copyobj(get(axesHandles, 'Children'), gca);
                title(get(axesHandles, 'Title').String);
                xlabel(get(axesHandles, 'XLabel').String);
                ylabel(get(axesHandles, 'YLabel').String);
                colorbar
                close(figHandles(i));
            end

        end

        function dehomogenize(obj)
            a = obj.harmonicVector;
            s.nCells             = 75;
            s.cellLevelSetParams = obj.createLevelSetCellParams();
            s.mesh               = obj.mesh;
            s.orientationA       = obj.harmonicVector;
            d = Dehomogenizer(s);
            d.compute();
            d.plot();
        end

        function s = createLevelSetCellParams(obj)
            s.type   = 'smoothRectangle';
            % s.type   = 'rectangleInclusion';
            s.xSide = obj.createFunction(obj.experimentData.dataRes.DesignVar1);%0.85*ones(size(obj.mesh.coord,1),1);
            s.ySide = obj.createFunction(obj.experimentData.dataRes.DesignVar2);%0.85*ones(size(obj.mesh.coord,1),1);
            s.pnorm  = obj.createFunction(obj.experimentData.dataRes.SuperEllipseExponent);
            s.ndim   = 2;
        end

        function f = createFunction(obj,value)
            s.fValues = value;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            f = LagrangianFunction(s);
            f = project(f,'P1D');
        end        



        %%%%%%%%%%%%%%%%%%%%%%%%


        function aBar = createOrientationByHand(obj)
            % alpha = obj.createLinealWithNoiseAngle();
            alpha = obj.createRadialWithNoiseAngle();
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