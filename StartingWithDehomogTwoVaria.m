% classdef StartingWithDehomogTwoVaria < handle
% 
%     properties (Access = private)
%         mesh
%         bField
%         rhoField
%     end
% 
%     methods (Access = public)
% 
%         function obj = StartingWithDehomogTwoVaria(bDesignVariable, rhoDesignVariable)
%             obj.init()
%             obj.mesh     = bDesignVariable.mesh;
%             obj.bField   = bDesignVariable;
%             obj.rhoField = rhoDesignVariable;
%             obj.createLevelSet(0.035)
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function init(obj)
%         end
% 
%         function ls = createLevelSet(obj, eps)
%             close all
% 
%             x = obj.mesh.coordElem;
%             fValues = obj.geometricalFunction(x, eps);
% 
%             s.mesh    = obj.mesh;
%             s.fValues = fValues(:);
%             s.order   = 'P1D';
%             ls = LagrangianFunction(s);
% 
%             sF.trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
%             sF.mesh  = obj.mesh;
%             filter   = FilterLump(sF);
%             ls       = filter.compute(ls, 3);
% 
%             sUm.backgroundMesh = obj.mesh;
%             sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
%             uMesh              = UnfittedMesh(sUm);
%             uMesh.compute(-ls.fValues);
%             uMesh.plot()
% 
% 
% 
% 
%         end
% 
%         % function fH = geometricalFunction(obj, xV, eps)
%         %     s.fHandle = @(x) obj.coordFun(x);
%         %     s.ndimf   = 2;
%         %     s.mesh    = obj.mesh;
%         %     xCoord    = AnalyticalFunction(s);
%         % 
%         %     s.operation = @(xV) obj.evaluateCellCoord(xV, eps, xCoord);
%         %     s.mesh      = obj.mesh;
%         %     s.ndimf     = xCoord.ndimf;
%         %     txi         = DomainFunction(s);
%         %     xiV         = txi.evaluate(xV);            
%         % 
%         % 
%         %     b   = obj.bField.evaluate(xV);
%         %     rho = obj.rhoField.evaluate(xV);
%         % 
%         % 
%         %     a = exp(b.^2);
%         %     d = (1 + b.^2) ./ a;
%         % 
%         % 
%         %     xi1  = xiV(1,:,:);
%         %     xi2  = xiV(2,:,:);
%         %     xic1 = xi1 - 0.5;
%         %     xic2 = xi2 - 0.5;
%         % 
%         % 
%         %     phi1 = a.*xic1 + b.*xic2;
%         %     phi2 = b.*xic1 + d.*xic2;
%         % 
%         %     % rhoMin = 0.1;
%         %     % rhoMax = 0.7;
%         %     % rhoEff = min(max(rho,rhoMin),rhoMax);
%         %
%         %     fH = max(abs(phi1)./(rho/2),abs(phi2)./(rho/2)) - 1;
%         % end
%         function fH = geometricalFunction(obj, xV, eps)
%             s.fHandle   = @(x) obj.coordFun(x);
%             s.ndimf     = 2;
%             s.mesh      = obj.mesh;
%             xCoord      = AnalyticalFunction(s);
% 
%             s.operation = @(xV) obj.evaluateCellCoord(xV, eps, xCoord);
%             s.mesh      = obj.mesh;
%             s.ndimf     = xCoord.ndimf;
%             txi         = DomainFunction(s);
%             xiV         = txi.evaluate(xV);
% 
%             b   = obj.bField.evaluate(xV);
%             rho = obj.rhoField.evaluate(xV);
% 
%             % Projecta rho para 0/1 com threshold 0.5
%             % Suaviza a transicao com sigmoid para evitar descontinuidades
%             beta     = 50;                          % controla a nitidez
%             rhoProj  = 1 ./ (1 + exp(-beta*(rho - 0.5)));
%             rhoProj  = min(max(real(rhoProj), 0), 1 - 1e-6);
% 
%             a = exp(b.^2);
%             d = (1 + b.^2) ./ a;
% 
%             xic1 = xiV(1,:,:) ;
%             xic2 = xiV(2,:,:) ;
% 
%             e1 =  d.*xic1 - b.*xic2;
%             e2 = -b.*xic1 + a.*xic2;
% 
%             holeSize = sqrt(1 - rhoProj) / 2;
%             holeSize = max(holeSize, 1e-6);
% 
%             fH = max(abs(e1)./holeSize, abs(e2)./holeSize) - 1;
%         end
% 
%         function f = coordFun(obj, x)
%             x1 = x(1,:,:);
%             x2 = x(2,:,:);
%             f  = [x1; x2];
%         end
% 
%         function txi = evaluateCellCoord(obj, xV, eps, xCoord)
%             x   = xCoord.evaluate(xV);
%             y   = obj.computeMicroCoordinate(x, eps);
%             txi = obj.periodicFunction(y);
%         end
% 
%         function y = computeMicroCoordinate(obj, x, eps)
%             y = x / eps(1);
%         end
% 
%     end
% 
%     methods (Access = private, Static)
% 
%         function f = periodicFunction(y)
%             f = y - floor(y);
%             f = f - 0.5;
%         end
% 
%     end
% 
% end
classdef StartingWithDehomogTwoVaria < handle

    properties (Access = private)

        mesh
        bField
        rhoField

        epsMicro
        NxPhase
        NyPhase

        phase1Interpolator
        phase2Interpolator
        rhoInterpolator

        phaseData
        levelSet

    end

    methods (Access = public)

        %% =========================================================
        % CONSTRUCTOR
        %
        % Important:
        % The constructor only stores the fields.
        % It does NOT solve phases, generate figures or create
        % UnfittedMesh automatically.
        %% =========================================================

        function obj = StartingWithDehomogTwoVaria( ...
                bDesignVariable, rhoDesignVariable, ...
                epsMicro, NxPhase, NyPhase)

            if nargin < 3 || isempty(epsMicro)
                epsMicro = 0.10;
            end

            if nargin < 4 || isempty(NxPhase)
                NxPhase = 200;
            end

            if nargin < 5 || isempty(NyPhase)
                NyPhase = 180;
            end

            obj.mesh     = bDesignVariable.mesh;
            obj.bField   = bDesignVariable;
            obj.rhoField = rhoDesignVariable;

            obj.epsMicro = epsMicro;
            obj.NxPhase  = NxPhase;
            obj.NyPhase  = NyPhase;

            obj.phaseData = struct();
            obj.levelSet  = [];

            fprintf('\n');
            fprintf('============================================================\n');
            fprintf('SAFE COMPATIBLE-PHASE DEHOMOGENIZATION\n');
            fprintf('============================================================\n');
            fprintf('epsMicro = %.6e\n',obj.epsMicro);
            fprintf('phase grid = %d x %d\n',obj.NxPhase,obj.NyPhase);
            fprintf('\nNo phase system has been solved yet.\n');
            fprintf('Next command:\n');
            fprintf('    dehomo.checkResolution()\n');
            fprintf('============================================================\n');

        end

        %% =========================================================
        % PUBLIC STEP 1: RESOLUTION CHECK
        %% =========================================================

        function resolution = checkResolution(obj)

            coord = obj.mesh.coord;

            xUnique = unique(coord(:,1));
            yUnique = unique(coord(:,2));

            dxValues = diff(xUnique);
            dyValues = diff(yUnique);

            dxValues = dxValues(dxValues > 1e-12);
            dyValues = dyValues(dyValues > 1e-12);

            if isempty(dxValues) || isempty(dyValues)
                error('Could not estimate the background mesh size.');
            end

            hMesh = min([dxValues(:);dyValues(:)]);

            rhoValues = real(obj.rhoField.fValues(:));

            % Only consider the region where material is actually present.
            rhoActive = rhoValues(rhoValues > 0.05);

            if isempty(rhoActive)
                rhoActive = rhoValues;
            end

            rhoMedian = median(rhoActive);
            rhoMean   = mean(rhoActive);

            fprintf('Active nodes used in resolution check = %d / %d\n', ...
                numel(rhoActive),numel(rhoValues));

            rhoMedian = min(max(rhoMedian,1e-6),1-1e-6);
            rhoMean   = min(max(rhoMean,1e-6),1-1e-6);

            relativeWallMedian = ...
                (1-sqrt(1-rhoMedian))/2;

            relativeWallMean = ...
                (1-sqrt(1-rhoMean))/2;

            wallMedian = obj.epsMicro*relativeWallMedian;
            wallMean   = obj.epsMicro*relativeWallMean;

            nElemMedian = wallMedian/hMesh;
            nElemMean   = wallMean/hMesh;

            resolution = struct();

            resolution.hMesh       = hMesh;
            resolution.rhoMedian   = rhoMedian;
            resolution.rhoMean     = rhoMean;
            resolution.wallMedian  = wallMedian;
            resolution.wallMean    = wallMean;
            resolution.nElemMedian = nElemMedian;
            resolution.nElemMean   = nElemMean;

            fprintf('\n');
            fprintf('============================================================\n');
            fprintf('GEOMETRIC RESOLUTION CHECK\n');
            fprintf('============================================================\n');
            fprintf('Background-mesh size h       = %.6e\n',hMesh);
            fprintf('Median rho                   = %.6e\n',rhoMedian);
            fprintf('Mean rho                     = %.6e\n',rhoMean);
            fprintf('Wall thickness at median rho = %.6e\n',wallMedian);
            fprintf('Wall thickness at mean rho   = %.6e\n',wallMean);
            fprintf('Elements/wall, median rho    = %.3f\n',nElemMedian);
            fprintf('Elements/wall, mean rho      = %.3f\n',nElemMean);

            if min(nElemMedian,nElemMean) < 2
                fprintf('\nSTATUS: SEVERELY UNDER-RESOLVED\n');
                fprintf('Increase epsMicro or use a finer reconstruction mesh.\n');
            elseif min(nElemMedian,nElemMean) < 3
                fprintf('\nSTATUS: MARGINALLY RESOLVED\n');
                fprintf('Three to four elements through the wall are recommended.\n');
            else
                fprintf('\nSTATUS: ACCEPTABLE FOR AN INITIAL TEST\n');
            end

            fprintf('============================================================\n');

        end

        %% =========================================================
        % PUBLIC STEP 2: COMPUTE PHASES ONLY
        %% =========================================================

        function computePhases(obj)

            fprintf('\n');
            fprintf('============================================================\n');
            fprintf('COMPUTING COMPATIBLE PHASE FIELDS\n');
            fprintf('============================================================\n');
            fprintf('Grid: %d x %d\n',obj.NxPhase,obj.NyPhase);

            coord = obj.mesh.coord;

            xMin = min(coord(:,1));
            xMax = max(coord(:,1));
            yMin = min(coord(:,2));
            yMax = max(coord(:,2));

            xGrid = linspace(xMin,xMax,obj.NxPhase);
            yGrid = linspace(yMin,yMax,obj.NyPhase);

            [X,Y] = meshgrid(xGrid,yGrid);

            dx = xGrid(2)-xGrid(1);
            dy = yGrid(2)-yGrid(1);

            %% ------------------------------------------------------------
            % Interpolation of the recovered P1 fields
            %% ------------------------------------------------------------

            fprintf('Interpolating nodal b and rho on the phase grid...\n');

            coordNodes = obj.mesh.coord;

            bNodal   = real(obj.bField.fValues(:));
            rhoNodal = real(obj.rhoField.fValues(:));

            if size(coordNodes,1) ~= numel(bNodal)
                error(['The recovered b field is not nodal P1: ', ...
                    'nodes = %d, b values = %d.'], ...
                    size(coordNodes,1),numel(bNodal));
            end

            if size(coordNodes,1) ~= numel(rhoNodal)
                error(['The recovered rho field is not nodal P1: ', ...
                    'nodes = %d, rho values = %d.'], ...
                    size(coordNodes,1),numel(rhoNodal));
            end

            Fb = scatteredInterpolant( ...
                coordNodes(:,1), ...
                coordNodes(:,2), ...
                bNodal, ...
                'linear', ...
                'nearest');

            Frho = scatteredInterpolant( ...
                coordNodes(:,1), ...
                coordNodes(:,2), ...
                rhoNodal, ...
                'linear', ...
                'nearest');

            bValues   = Fb(X,Y);
            rhoValues = Frho(X,Y);

            bValues   = real(bValues);
            rhoValues = real(rhoValues);

            bValues = min(max(bValues,-0.8),0.8);
            rhoValues = min(max(rhoValues,1e-6),1-1e-6);

            fprintf('Interpolated values:\n');
            fprintf('  b   min/max = %.6e %.6e\n', ...
                min(bValues(:)),max(bValues(:)));
            fprintf('  rho min/max = %.6e %.6e\n', ...
                min(rhoValues(:)),max(rhoValues(:)));

            % Light grid-scale smoothing of b.
            bValues = obj.smoothScalarGrid(bValues,1);

            %% ------------------------------------------------------------
            % Reciprocal lattice fields
            %% ------------------------------------------------------------

            aValues = exp(bValues.^2);
            dValues = (1+bValues.^2)./aValues;

            q1x =  dValues/obj.epsMicro;
            q1y = -bValues/obj.epsMicro;

            q2x = -bValues/obj.epsMicro;
            q2y =  aValues/obj.epsMicro;

            fprintf('Solving phase 1 with LSQR...\n');

            phi1 = obj.solveGradientProjectionLSQR( ...
                q1x,q1y,dx,dy);

            fprintf('Solving phase 2 with LSQR...\n');

            phi2 = obj.solveGradientProjectionLSQR( ...
                q2x,q2y,dx,dy);

            %% ------------------------------------------------------------
            % Compatibility diagnostics
            %% ------------------------------------------------------------

            fprintf('Computing compatibility diagnostics...\n');

            [dq1x_dy,~] = gradient(q1x,dy,dx);
            [~,dq1y_dx] = gradient(q1y,dy,dx);

            [dq2x_dy,~] = gradient(q2x,dy,dx);
            [~,dq2y_dx] = gradient(q2y,dy,dx);

            curl1 = dq1y_dx-dq1x_dy;
            curl2 = dq2y_dx-dq2x_dy;

            [phi1Y,phi1X] = gradient(phi1,dy,dx);
            [phi2Y,phi2X] = gradient(phi2,dy,dx);

            residual1 = sqrt( ...
                (phi1X-q1x).^2 + ...
                (phi1Y-q1y).^2);

            residual2 = sqrt( ...
                (phi2X-q2x).^2 + ...
                (phi2Y-q2y).^2);

            curl1RMS = sqrt(mean(curl1(:).^2));
            curl2RMS = sqrt(mean(curl2(:).^2));

            residual1RMS = sqrt(mean(residual1(:).^2));
            residual2RMS = sqrt(mean(residual2(:).^2));

            fprintf('\n');
            fprintf('RMS curl 1     = %.6e\n',curl1RMS);
            fprintf('RMS curl 2     = %.6e\n',curl2RMS);
            fprintf('RMS residual 1 = %.6e\n',residual1RMS);
            fprintf('RMS residual 2 = %.6e\n',residual2RMS);

            %% ------------------------------------------------------------
            % Interpolators for reconstruction
            %% ------------------------------------------------------------

            obj.phase1Interpolator = ...
                griddedInterpolant( ...
                {yGrid,xGrid},phi1, ...
                'linear','nearest');

            obj.phase2Interpolator = ...
                griddedInterpolant( ...
                {yGrid,xGrid},phi2, ...
                'linear','nearest');

            obj.rhoInterpolator = ...
                griddedInterpolant( ...
                {yGrid,xGrid},rhoValues, ...
                'linear','nearest');

            obj.phaseData.X         = X;
            obj.phaseData.Y         = Y;
            obj.phaseData.b         = bValues;
            obj.phaseData.rho       = rhoValues;
            obj.phaseData.phi1      = phi1;
            obj.phaseData.phi2      = phi2;
            obj.phaseData.curl1     = curl1;
            obj.phaseData.curl2     = curl2;
            obj.phaseData.residual1 = residual1;
            obj.phaseData.residual2 = residual2;

            obj.phaseData.curl1RMS     = curl1RMS;
            obj.phaseData.curl2RMS     = curl2RMS;
            obj.phaseData.residual1RMS = residual1RMS;
            obj.phaseData.residual2RMS = residual2RMS;

            fprintf('\nCompatible phases computed successfully.\n');
            fprintf('No level set or UnfittedMesh was created.\n');
            fprintf('Next command:\n');
            fprintf('    dehomo.savePhaseResults(''CompatiblePhaseResults_60x30.mat'')\n');
            fprintf('============================================================\n');

        end

        %% =========================================================
        % PUBLIC STEP 3: PLOT ONLY ONE DIAGNOSTIC
        %% =========================================================

        function plotOneDiagnostic(obj,name)

            if isempty(fieldnames(obj.phaseData))
                error('Run dehomo.computePhases() first.');
            end

            X = obj.phaseData.X;
            Y = obj.phaseData.Y;

            switch lower(name)

                case 'b'
                    Z = obj.phaseData.b;
                    plotTitle = 'Field $b(x)$';

                case 'rho'
                    Z = obj.phaseData.rho;
                    plotTitle = 'Field $\rho(x)$';

                case 'phi1'
                    Z = obj.phaseData.phi1;
                    plotTitle = 'Compatible phase $\phi_1(x)$';

                case 'phi2'
                    Z = obj.phaseData.phi2;
                    plotTitle = 'Compatible phase $\phi_2(x)$';

                case 'curl1'
                    Z = obj.phaseData.curl1;
                    plotTitle = ...
                        'Compatibility defect $\mathrm{curl}(q_1)$';

                case 'curl2'
                    Z = obj.phaseData.curl2;
                    plotTitle = ...
                        'Compatibility defect $\mathrm{curl}(q_2)$';

                case 'residual1'
                    Z = obj.phaseData.residual1;
                    plotTitle = 'Phase-projection residual 1';

                case 'residual2'
                    Z = obj.phaseData.residual2;
                    plotTitle = 'Phase-projection residual 2';

                otherwise
                    error(['Unknown diagnostic "%s". Use b, rho, ', ...
                        'phi1, phi2, curl1, curl2, residual1 or residual2.'], ...
                        name);
            end

            figure('Color','w');
            contourf(X,Y,Z,30,'LineColor','none');
            axis equal tight;
            colorbar;
            xlabel('$x$','Interpreter','latex');
            ylabel('$y$','Interpreter','latex');
            title(plotTitle,'Interpreter','latex');

        end

        %% =========================================================
        % PUBLIC STEP 4: SAVE PHASE DATA BEFORE LEVEL-SET CREATION
        %% =========================================================

        function savePhaseResults(obj,fileName)

            if nargin < 2 || isempty(fileName)
                fileName = 'CompatiblePhaseResults.mat';
            end

            if isempty(fieldnames(obj.phaseData))
                error('Run dehomo.computePhases() first.');
            end

            phaseResults = obj.phaseData;
            epsMicro     = obj.epsMicro;
            NxPhase      = obj.NxPhase;
            NyPhase      = obj.NyPhase;

            save(fileName, ...
                'phaseResults','epsMicro','NxPhase','NyPhase','-v7.3');

            fprintf('Phase results saved in:\n%s\n',fileName);

        end

        %% =========================================================
        % PUBLIC STEP 5: BUILD LEVEL-SET VALUES ONLY
        %
        % Does not create UnfittedMesh yet.
        %% =========================================================

        function ls = buildLevelSet(obj)

            if isempty(obj.phase1Interpolator) || ...
                    isempty(obj.phase2Interpolator)

                error('Run dehomo.computePhases() first.');
            end

            fprintf('\nBuilding level-set values...\n');

            x = obj.mesh.coordElem;
            fValues = obj.evaluateGeometricalFunction(x);

            s.mesh    = obj.mesh;
            s.fValues = fValues(:);
            s.order   = 'P1D';

            ls = LagrangianFunction(s);

            obj.levelSet = ls;

            fprintf('Level-set values generated successfully.\n');
            fprintf('No UnfittedMesh was created.\n');

        end

        %% =========================================================
        % PUBLIC STEP 6: PLOT LEVEL-SET FUNCTION
        %% =========================================================

        function plotLevelSet(obj)

            if isempty(obj.levelSet)
                error('Run dehomo.buildLevelSet() first.');
            end

            figure('Color','w');
            obj.levelSet.plot();

            title(sprintf( ...
                'Compatible-phase level set, epsilon = %.4f', ...
                obj.epsMicro), ...
                'Interpreter','none');

        end

        %% =========================================================
        % PUBLIC STEP 7: CREATE UNFITTED MESH
        %
        % Run only after phase and level-set files have been saved.
        %% =========================================================

        function createAndPlotUnfittedMesh(obj,useNegativeSign)

            if nargin < 2
                useNegativeSign = false;
            end

            if isempty(obj.levelSet)
                error('Run dehomo.buildLevelSet() first.');
            end

            fprintf('\nCreating UnfittedMesh...\n');

            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;

            uMesh = UnfittedMesh(sUm);

            if useNegativeSign
                uMesh.compute(-obj.levelSet.fValues);
            else
                uMesh.compute(obj.levelSet.fValues);
            end

            uMesh.plot();

            title(sprintf( ...
                'Compatible-phase dehomogenization, epsilon = %.4f', ...
                obj.epsMicro), ...
                'Interpreter','none');

        end

        %% =========================================================
        % PUBLIC SUMMARY
        %% =========================================================

        function printSummary(obj)

            fprintf('\n');
            fprintf('============================================================\n');
            fprintf('DEHOMOGENIZATION OBJECT SUMMARY\n');
            fprintf('============================================================\n');
            fprintf('epsMicro  = %.6e\n',obj.epsMicro);
            fprintf('NxPhase   = %d\n',obj.NxPhase);
            fprintf('NyPhase   = %d\n',obj.NyPhase);
            fprintf('Phases computed: %d\n', ...
                ~isempty(obj.phase1Interpolator));
            fprintf('Level set built: %d\n', ...
                ~isempty(obj.levelSet));

            if ~isempty(fieldnames(obj.phaseData))
                fprintf('RMS curl 1     = %.6e\n', ...
                    obj.phaseData.curl1RMS);
                fprintf('RMS curl 2     = %.6e\n', ...
                    obj.phaseData.curl2RMS);
                fprintf('RMS residual 1 = %.6e\n', ...
                    obj.phaseData.residual1RMS);
                fprintf('RMS residual 2 = %.6e\n', ...
                    obj.phaseData.residual2RMS);
            end

            fprintf('============================================================\n');

        end

    end

    methods (Access = private)

        %% =========================================================
        % MEMORY-SAFE LSQR PROJECTION
        %
        % min_phi ||G phi - q||_2
        %% =========================================================

        function phi = solveGradientProjectionLSQR( ...
                obj,qx,qy,dx,dy) %#ok<INUSL>

            [Ny,Nx] = size(qx);

            nNodes = Nx*Ny;

            nHorizontalEdges = Ny*(Nx-1);
            nVerticalEdges   = (Ny-1)*Nx;

            nRows = ...
                nHorizontalEdges + ...
                nVerticalEdges + 1;

            maxNonZeros = ...
                2*(nHorizontalEdges+nVerticalEdges)+1;

            rowIndex = zeros(maxNonZeros,1);
            colIndex = zeros(maxNonZeros,1);
            values   = zeros(maxNonZeros,1);

            rhs = zeros(nRows,1);

            row = 0;
            nz  = 0;

            % Horizontal edges
            for iy = 1:Ny
                for ix = 1:Nx-1

                    row = row+1;

                    nodeLeft = sub2ind([Ny,Nx],iy,ix);
                    nodeRight = sub2ind([Ny,Nx],iy,ix+1);

                    nz = nz+1;
                    rowIndex(nz) = row;
                    colIndex(nz) = nodeLeft;
                    values(nz)   = -1/dx;

                    nz = nz+1;
                    rowIndex(nz) = row;
                    colIndex(nz) = nodeRight;
                    values(nz)   = 1/dx;

                    rhs(row) = ...
                        0.5*(qx(iy,ix)+qx(iy,ix+1));

                end
            end

            % Vertical edges
            for iy = 1:Ny-1
                for ix = 1:Nx

                    row = row+1;

                    nodeBottom = sub2ind([Ny,Nx],iy,ix);
                    nodeTop = sub2ind([Ny,Nx],iy+1,ix);

                    nz = nz+1;
                    rowIndex(nz) = row;
                    colIndex(nz) = nodeBottom;
                    values(nz)   = -1/dy;

                    nz = nz+1;
                    rowIndex(nz) = row;
                    colIndex(nz) = nodeTop;
                    values(nz)   = 1/dy;

                    rhs(row) = ...
                        0.5*(qy(iy,ix)+qy(iy+1,ix));

                end
            end

            % Anchor to remove the arbitrary additive constant
            anchorWeight = 1e3;

            row = row+1;
            nz  = nz+1;

            rowIndex(nz) = row;
            colIndex(nz) = 1;
            values(nz)   = anchorWeight;

            rhs(row) = 0;

            rowIndex = rowIndex(1:nz);
            colIndex = colIndex(1:nz);
            values   = values(1:nz);

            G = sparse( ...
                rowIndex,colIndex,values,nRows,nNodes);

            fprintf('System dimensions: %d x %d\n', ...
                size(G,1),size(G,2));

            fprintf('Nonzero entries: %d\n',nnz(G));

            tolLSQR = 1e-6;
            maxIter = 300;

            [phiVector,flag,relres,iter] = lsqr( ...
                G,rhs,tolLSQR,maxIter);

            fprintf('LSQR flag       = %d\n',flag);
            fprintf('LSQR relres     = %.6e\n',relres);
            fprintf('LSQR iterations = %d\n',iter);

            if flag ~= 0 && flag ~= 1
                warning('LSQR returned flag %d.',flag);
            end

            phi = reshape(phiVector,Ny,Nx);

            % Explicitly release temporary sparse data.
            clear G rhs rowIndex colIndex values phiVector

        end

        %% =========================================================
        % LEVEL-SET GEOMETRY
        %% =========================================================

        function fH = evaluateGeometricalFunction(obj,xV)

            s.fHandle = @(x) obj.coordFun(x);
            s.ndimf   = 2;
            s.mesh    = obj.mesh;

            xCoord = AnalyticalFunction(s);
            x = xCoord.evaluate(xV);

            x1 = x(1,:,:);
            x2 = x(2,:,:);

            phi1 = obj.phase1Interpolator(x2,x1);
            phi2 = obj.phase2Interpolator(x2,x1);

            xi1 = phi1-floor(phi1)-0.5;
            xi2 = phi2-floor(phi2)-0.5;

            rho = obj.rhoInterpolator(x2,x1);

            rho = min(max(real(rho),1e-6),1-1e-6);

            holeHalfSize = sqrt(1-rho)/2;
            holeHalfSize = max(holeHalfSize,1e-6);

            fH = max( ...
                abs(xi1)./holeHalfSize, ...
                abs(xi2)./holeHalfSize)-1;

        end

        %% =========================================================
        % LIGHT GRID SMOOTHING
        %% =========================================================

        function uSmooth = smoothScalarGrid(obj,u,nSteps) %#ok<INUSL>

            uSmooth = u;

            kernel = [
                0 1 0;
                1 4 1;
                0 1 0
                ];

            kernel = kernel/sum(kernel(:));

            for k = 1:nSteps

                oldValues = uSmooth;

                uSmooth = conv2( ...
                    oldValues,kernel,'same');

                uSmooth(1,:)   = oldValues(1,:);
                uSmooth(end,:) = oldValues(end,:);
                uSmooth(:,1)   = oldValues(:,1);
                uSmooth(:,end) = oldValues(:,end);

            end

        end

        %% =========================================================
        % COORDINATES
        %% =========================================================

        function f = coordFun(obj,x) %#ok<INUSL>

            x1 = x(1,:,:);
            x2 = x(2,:,:);

            f = [x1;x2];

        end

    end

end