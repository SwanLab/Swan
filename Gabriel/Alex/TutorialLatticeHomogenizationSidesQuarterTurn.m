classdef TutorialLatticeHomogenizationSidesQuarterTurn < handle
    properties (Access = public)
        dataset
    end
    properties (Access = private)
        E
        nu
        divUnit
        mU
        cellTransformation
        wallAngleMin
        wallAngleMax
        nWallAngles
        nSides
        nWorkers
        folder
    end
    properties (Access = private)
        settings
        generator
    end
    methods (Access = public)
        function obj = TutorialLatticeHomogenizationSidesQuarterTurn()
            obj.init();
            obj.createSettings();
            obj.createGenerator();
        end
        function compute(obj)
            obj.dataset = generateSidesLatticeDataset(obj.settings);
        end
        function plotMicrostructures(obj,aValues,m1,m2,nCells,resolution)
            if nargin < 5
                nCells = 3;
            end
            if nargin < 6
                resolution = 120;
            end
            figure('Color','w')
            tiledlayout(1,numel(aValues),'TileSpacing','compact','Padding','compact')
            for i = 1:numel(aValues)
                [X,Y,isSolid,T] = obj.generator.sampleMicrostructure(aValues(i),m1,m2,nCells,resolution);
                nexttile
                pcolor(X,Y,double(isSolid))
                shading flat
                colormap(gca,[1 1 1; 0.15 0.15 0.15])
                clim([0 1])
                axis equal off
                title({sprintf('a = %.3g, m_1 = %.2g, m_2 = %.2g, \\rho = %.2g',aValues(i),m1,m2,1-m1*m2), ...
                       sprintf('wall at %.1f^\\circ, det T = %.4f',atan2d(T(2,1),T(1,1)),det(T))},'FontWeight','normal')
            end
        end
        function plotTensor(obj,m1,m2)
            d           = obj.loadDataset();
            isComponent = ~ismember(d.columns,["m1","m2","a","rho"]);
            sides       = d.data(:,ismember(d.columns,["m1","m2"]));
            [~,i]       = min(vecnorm(sides-[m1,m2],2,2));
            rows        = all(sides == sides(i,:),2);
            a           = d.data(rows,d.columns == "a");
            C           = d.data(rows,isComponent);
            names       = d.columns(isComponent);
            figure('Color','w')
            for k = 1:numel(names)
                subplot(2,3,k)
                plot(a,C(:,k),'.-','LineWidth',1.2)
                grid on
                xlabel('a')
                title(names(k))
            end
            sgtitle(sprintf('Homogenized tensor, m_1 = %.3g, m_2 = %.3g',sides(i,1),sides(i,2)))
        end
        function a = getAValues(obj)
            a = obj.settings.dataset.a;
        end
        function m = getSideValues(obj)
            m = obj.settings.dataset.m;
        end
    end
    methods (Access = private)
        function init(obj)
            obj.E                  = 1;
            obj.nu                 = 1/3;
            obj.divUnit            = 160;
            obj.mU                 = 0.99;
            obj.cellTransformation = 'QuarterTurn';
            obj.wallAngleMin       = -33.75;
            obj.wallAngleMax       = 123.75;
            obj.nWallAngles        = 29;
            obj.nSides             = 21;
            obj.nWorkers           = 4;
            obj.folder             = 'latticeSidesQuarterTurnDataset';
        end
        function createSettings(obj)
            s                    = latticeSidesSettings();
            s.material.E         = obj.E;
            s.material.nu        = obj.nu;
            s.mesh.divUnit       = obj.divUnit;
            s.hole.mU            = obj.mU;
            s.cellTransformation = obj.cellTransformation;
            s.dataset.m          = obj.computeSideValues();
            s.dataset.a          = obj.computeAValues();
            s.dataset.symmetries = {};
            s.dataset.nWorkers   = obj.nWorkers;
            s.dataset.folder     = obj.folder;
            obj.settings         = s;
        end
        function createGenerator(obj)
            obj.generator = SidesLatticeDatasetGenerator(obj.settings);
        end
        function m = computeSideValues(obj)
            m = obj.mU*(1-cos(linspace(0,pi,obj.nSides)))/2;
        end
        function a = computeAValues(obj)
            theta = linspace(obj.wallAngleMin,obj.wallAngleMax,obj.nWallAngles);
            a     = sind(theta)./(sind(theta)+cosd(theta));
        end
        function d = loadDataset(obj)
            if isempty(obj.dataset)
                obj.dataset = load(fullfile(obj.folder,obj.settings.dataset.outputFile));
            end
            d = obj.dataset;
        end
    end
end