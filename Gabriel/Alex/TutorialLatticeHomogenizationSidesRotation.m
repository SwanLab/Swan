classdef TutorialLatticeHomogenizationSidesRotation < handle
    properties (Access = public)
        dataset
    end
    properties (Access = private)
        E
        nu
        divUnit
        mU
        cellTransformation
        quarterTurnFile
        period
        nWorkers
        folder
        fullTurnFile
    end
    properties (Access = private)
        settings
        generator
    end
    methods (Access = public)
        function obj = TutorialLatticeHomogenizationSidesRotation()
            obj.init();
            obj.createSettings();
            obj.createGenerator();
        end
        function compute(obj)
            generateSidesLatticeDataset(obj.settings);
            obj.merge();
            obj.check();
        end
        function merge(obj)
            old     = obj.loadQuarterTurnDataset();
            new     = load(fullfile(obj.folder,obj.settings.dataset.outputFile));
            isAngle = old.columns == "a";
            if ~isequal(old.columns,new.columns)
                error('TutorialLatticeHomogenizationSidesRotation:DifferentColumns','The two datasets do not have the same columns.');
            end
            data            = [old.data;new.data];
            tensors         = cat(5,old.tensors,new.tensors);
            isParameter     = ismember(old.columns,["m1","m2","a"]);
            [~,order]       = sortrows(round(data(:,isParameter),10));
            full.data       = data(order,:);
            full.tensors    = tensors(:,:,:,:,order);
            full.columns    = old.columns;
            full.columns(isAngle) = "theta";
            full.period     = obj.period;
            full.settings   = obj.settings;
            obj.dataset     = full;
            save(obj.fullTurnFile,'-struct','full');
            theta = unique(round(full.data(:,isAngle),10));
            fprintf('Dataset with %d samples and %d angles from %.4g to %.4g degrees saved in %s\n',size(full.data,1),numel(theta),min(theta),max(theta),obj.fullTurnFile);
        end
        function report = check(obj,rhoMin)
            if nargin < 2
                rhoMin = 0.05;
            end
            d       = obj.loadDataset();
            m1      = round(d.data(:,d.columns == "m1"),10);
            m2      = round(d.data(:,d.columns == "m2"),10);
            theta   = round(d.data(:,d.columns == "theta"),10);
            rho     = d.data(:,d.columns == "rho");
            partner = round(theta-obj.period/2,10);
            [hasPartner,rows] = ismember([m2,m1,partner],[m1,m2,theta],'rows');
            isNew   = ismember(theta,round(obj.settings.dataset.a(:),10));
            keep    = hasPartner & rho > rhoMin;
            report.all = obj.compare('(m1,m2,theta) against (m2,m1,theta - 90), all angles with a partner',d.tensors,rows,keep,rho,theta);
            report.new = obj.compare('(m1,m2,theta) against (m2,m1,theta - 90), new angles only',d.tensors,rows,keep & isNew,rho,theta);
        end
        function plotMicrostructures(obj,thetaValues,m1,m2,nCells,resolution)
            if nargin < 5
                nCells = 3;
            end
            if nargin < 6
                resolution = 120;
            end
            figure('Color','w')
            tiledlayout(1,numel(thetaValues),'TileSpacing','compact','Padding','compact')
            for i = 1:numel(thetaValues)
                [X,Y,isSolid,T] = obj.generator.sampleMicrostructure(thetaValues(i),m1,m2,nCells,resolution);
                nexttile
                pcolor(X,Y,double(isSolid))
                shading flat
                colormap(gca,[1 1 1; 0.15 0.15 0.15])
                clim([0 1])
                axis equal off
                title({sprintf('\\theta = %.4g^\\circ, m_1 = %.2g, m_2 = %.2g, \\rho = %.2g',thetaValues(i),m1,m2,1-m1*m2), ...
                       sprintf('det T = %.4f',det(T))},'FontWeight','normal')
            end
        end
        function plotTensor(obj,m1,m2)
            d           = obj.loadDataset();
            isComponent = ~ismember(d.columns,["m1","m2","theta","rho"]);
            sides       = d.data(:,ismember(d.columns,["m1","m2"]));
            [~,i]       = min(vecnorm(sides-[m1,m2],2,2));
            rows        = all(sides == sides(i,:),2);
            theta       = d.data(rows,d.columns == "theta");
            C           = d.data(rows,isComponent);
            names       = d.columns(isComponent);
            figure('Color','w')
            for k = 1:numel(names)
                subplot(2,3,k)
                plot(theta,C(:,k),'.-','LineWidth',1.2)
                grid on
                xlabel('\theta (degrees)')
                title(names(k))
            end
            sgtitle(sprintf('Homogenized tensor, m_1 = %.3g, m_2 = %.3g',sides(i,1),sides(i,2)))
        end
        function theta = getNewAngles(obj)
            theta = obj.settings.dataset.a;
        end
        function file = getFullTurnFile(obj)
            file = obj.fullTurnFile;
        end
    end
    methods (Access = private)
        function init(obj)
            obj.E                  = 1;
            obj.nu                 = 1/3;
            obj.divUnit            = 160;
            obj.mU                 = 0.99;
            obj.cellTransformation = 'Rotation';
            obj.quarterTurnFile    = fullfile('latticeSidesQuarterTurnDataset','dataset.mat');
            obj.period             = 180;
            obj.nWorkers           = 4;
            obj.folder             = 'latticeSidesRotationDataset';
            obj.fullTurnFile       = fullfile(obj.folder,'datasetFullTurn.mat');
        end
        function createSettings(obj)
            old                  = obj.loadQuarterTurnDataset();
            s                    = latticeSidesSettings();
            s.material.E         = obj.E;
            s.material.nu        = obj.nu;
            s.mesh.divUnit       = obj.divUnit;
            s.hole.mU            = obj.mU;
            s.cellTransformation = obj.cellTransformation;
            s.dataset.m          = unique(round(old.data(:,old.columns == "m1"),10))';
            s.dataset.a          = obj.computeMissingAngles(old);
            s.dataset.symmetries = {};
            s.dataset.nWorkers   = obj.nWorkers;
            s.dataset.folder     = obj.folder;
            obj.settings         = s;
        end
        function createGenerator(obj)
            obj.generator = SidesLatticeDatasetGenerator(obj.settings);
        end
        function theta = computeMissingAngles(obj,old)
            existing = unique(round(old.data(:,old.columns == "a"),10));
            step     = existing(2) - existing(1);
            nMissing = round((min(existing) + obj.period - max(existing))/step) - 1;
            theta    = max(existing) + step*(1:nMissing);
        end
        function old = loadQuarterTurnDataset(obj)
            old     = load(obj.quarterTurnFile);
            isAngle = old.columns == "a";
            a       = old.data(:,isAngle);
            old.data(:,isAngle) = round(atan2d(a,1-a),10);
        end
        function d = loadDataset(obj)
            if isempty(obj.dataset)
                obj.dataset = load(obj.fullTurnFile);
            end
            d = obj.dataset;
        end
    end
    methods (Static, Access = private)
        function result = compare(name,C,rows,keep,rho,theta)
            index     = find(keep);
            C1        = C(:,:,:,:,index);
            C2        = C(:,:,:,:,rows(index));
            perSample = squeeze(sqrt(sum((C1-C2).^2,1:4)./sum(C2.^2,1:4)));
            [worst,i] = max(perSample);
            result.global = sqrt(sum((C1(:)-C2(:)).^2)/sum(C2(:).^2));
            result.max    = worst;
            fprintf('%s\n   %d samples | global %.3e | worst sample %.3e at rho %.3f, theta %.4g\n',name,numel(index),result.global,worst,rho(index(i)),theta(index(i)));
        end
    end
end