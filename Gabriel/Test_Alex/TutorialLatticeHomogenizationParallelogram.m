classdef TutorialLatticeHomogenizationParallelogram < handle
    properties (Access = public)
        dataset
    end
    properties (Access = private)
        E
        nu
        divUnit
        mU
        cellTransformation
        nAngles
        period
        sMax
        nBoosts
        rMax
        nStretches
        nSides
        nWorkers
        secondsPerGeometry
        folder
    end
    properties (Access = private)
        settings
        generator
        shapes
        isComputedShape
    end
    methods (Access = public)
        function obj = TutorialLatticeHomogenizationParallelogram()
            obj.init();
            obj.createShapes();
            obj.createSettings();
            obj.createGenerator();
        end
        function compute(obj)
            obj.dataset = generateCellLatticeDataset(obj.settings);
            obj.check();
        end
        function printSummary(obj)
            ds        = obj.settings.dataset;
            nGeom     = size(ds.cells,1)*numel(ds.m);
            computed  = obj.shapes(obj.isComputedShape,:);
            divisions = zeros(size(computed,1),2);
            for i = 1:size(computed,1)
                divisions(i,:) = obj.generator.computeDivisions([ds.a(1),computed(i,:)]);
            end
            disp(table(computed(:,1),computed(:,2),divisions(:,1),divisions(:,2),prod(divisions,2)/obj.divUnit^2, ...
                'VariableNames',{'s','r','divisions1','divisions2','meshSizeRelativeToSquare'}))
            fprintf('%d values of a from %.3f to %.3f, %d shapes (s,r) computed of %d, %d values of m1 and of m2\n', ...
                numel(ds.a),min(ds.a),max(ds.a),size(computed,1),size(obj.shapes,1),numel(ds.m));
            fprintf('%d geometries, %d samples computed, %d samples in the dataset\n', ...
                nGeom,nGeom*numel(ds.m),numel(ds.a)*size(obj.shapes,1)*numel(ds.m)^2);
            fprintf('Estimated time with %d workers at %.0f s per geometry: %.1f h\n', ...
                obj.nWorkers,obj.secondsPerGeometry,nGeom*obj.secondsPerGeometry/obj.nWorkers/3600);
        end
        function plotShapes(obj,a,m1,m2,nCells,resolution)
            if nargin < 5
                nCells = 2;
            end
            if nargin < 6
                resolution = 80;
            end
            sValues = unique(obj.shapes(:,1));
            rValues = unique(obj.shapes(:,2));
            figure('Color','w')
            tiledlayout(numel(rValues),numel(sValues),'TileSpacing','compact','Padding','compact')
            for j = numel(rValues):-1:1
                for i = 1:numel(sValues)
                    shape      = [sValues(i),rValues(j)];
                    isComputed = obj.isComputedShape(ismember(round(obj.shapes,10),round(shape,10),'rows'));
                    [X,Y,isSolid,T] = obj.generator.sampleMicrostructure([a,shape],m1,m2,nCells,resolution);
                    nexttile
                    pcolor(X,Y,double(isSolid))
                    shading flat
                    colormap(gca,[1 1 1; 0.15 0.15 0.15])
                    clim([0 1])
                    axis equal off
                    labels = {'twin','computed'};
                    title({sprintf('s = %.2g, r = %.3g (%s)',shape(1),shape(2),labels{isComputed+1}), ...
                           sprintf('det T = %.4f',det(T))},'FontWeight','normal')
                end
            end
            sgtitle(sprintf('a = %.3g, m_1 = %.2g, m_2 = %.2g',a,m1,m2))
        end
        function report = check(obj,rhoMin)
            if nargin < 2
                rhoMin = 0.05;
            end
            d       = obj.loadDataset();
            params  = d.data(:,ismember(d.columns,["m1","m2","a","s","r"]));
            rho     = d.data(:,d.columns == "rho");
            twin    = TurnBoostStretchTransformation().computeTwin(params);
            [hasTwin,rows] = ismember(round(twin,6),round(params,6),'rows');
            keep    = hasTwin & d.source == "computed" & d.source(max(rows,1)) == "computed" & rho > rhoMin;
            index   = find(keep);
            C1      = d.tensors(:,:,:,:,index);
            C2      = d.tensors(:,:,:,:,rows(index));
            perSample     = squeeze(sqrt(sum((C1-C2).^2,1:4)./sum(C2.^2,1:4)));
            [worst,i]     = max(perSample);
            report.global = sqrt(sum((C1(:)-C2(:)).^2)/sum(C2(:).^2));
            report.max    = worst;
            fprintf('Cells computed twice, (m1,m2,a,s,r) against (m2,m1,a turned 90 degrees,-s,1/r)\n');
            fprintf('   %d samples | global %.3e | worst sample %.3e at rho %.3f, a %.3f\n',numel(index),report.global,worst,rho(index(i)),params(index(i),3));
            eigMin = inf;
            for k = 1:size(d.tensors,5)
                C      = d.tensors(:,:,:,:,k);
                M      = [C(1,1,1,1),C(1,1,2,2),sqrt(2)*C(1,1,1,2);C(1,1,2,2),C(2,2,2,2),sqrt(2)*C(2,2,1,2);sqrt(2)*C(1,1,1,2),sqrt(2)*C(2,2,1,2),2*C(1,2,1,2)];
                eigMin = min(eigMin,min(eig((M+M')/2)));
            end
            report.minEigenvalue = eigMin;
            fprintf('   density against 1 - m1*m2: largest difference %.3e\n',max(abs(rho-(1-params(:,1).*params(:,2)))));
            fprintf('   smallest eigenvalue of the tensor over the dataset: %.3e\n',eigMin);
        end
        function a = getAValues(obj)
            a = obj.settings.dataset.a;
        end
        function shapes = getShapes(obj)
            shapes = obj.shapes;
        end
    end
    methods (Access = private)
        function init(obj)
            obj.E                  = 1;
            obj.nu                 = 1/3;
            obj.divUnit            = 160;
            obj.mU                 = 0.99;
            obj.cellTransformation = 'TurnBoostStretch';
            obj.nAngles            = 16;
            obj.period             = 180;
            obj.sMax               = 0.6;
            obj.nBoosts            = 5;
            obj.rMax               = 1.5;
            obj.nStretches         = 5;
            obj.nSides             = 21;
            obj.nWorkers           = 4;
            obj.secondsPerGeometry = 98;
            obj.folder             = 'latticeParallelogramDataset';
        end
        function createShapes(obj)
            sValues     = linspace(-obj.sMax,obj.sMax,obj.nBoosts);
            rValues     = obj.rMax.^linspace(-1,1,obj.nStretches);
            [S,R]       = ndgrid(sValues,rValues);
            obj.shapes  = [S(:),R(:)];
            tol         = 1e-10;
            obj.isComputedShape = S(:) > tol | (abs(S(:)) <= tol & R(:) >= 1 - tol);
        end
        function createSettings(obj)
            step                 = obj.period/obj.nAngles;
            wallAngles           = -obj.period/4 + step/2 + step*(0:obj.nAngles-1);
            computed             = obj.shapes(obj.isComputedShape,:);
            [iA,iShape]          = ndgrid(1:obj.nAngles,1:size(computed,1));
            s                    = latticeSidesSettings();
            s.material.E         = obj.E;
            s.material.nu        = obj.nu;
            s.mesh.divUnit       = obj.divUnit;
            s.hole.mU            = obj.mU;
            s.cellTransformation = obj.cellTransformation;
            s.dataset.m          = obj.mU*(1-cos(linspace(0,pi,obj.nSides)))/2;
            s.dataset.a          = sind(wallAngles)./(sind(wallAngles)+cosd(wallAngles));
            s.dataset.cells      = [s.dataset.a(iA(:))',computed(iShape(:),:)];
            s.dataset.symmetries = {};
            s.dataset.nWorkers   = obj.nWorkers;
            s.dataset.folder     = obj.folder;
            obj.settings         = s;
        end
        function createGenerator(obj)
            obj.generator = SidesLatticeDatasetGenerator(obj.settings);
        end
        function d = loadDataset(obj)
            if isempty(obj.dataset)
                obj.dataset = load(fullfile(obj.folder,obj.settings.dataset.outputFile));
            end
            d = obj.dataset;
        end
    end
end