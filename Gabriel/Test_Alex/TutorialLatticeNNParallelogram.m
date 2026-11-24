classdef TutorialLatticeNNParallelogram < handle
    properties (Access = public)
        surrogate
        report
    end
    properties (Access = private)
        datasetFile
        testFraction
        seed
        hiddenLayers
        learningRate
        maxEpochs
        earlyStop
        prefix
        outputFile
    end
    properties (Access = private)
        settings
        featureFile
        featureNames
        parameterNames
    end
    methods (Access = public)
        function obj = TutorialLatticeNNParallelogram()
            obj.init();
            obj.createFeatureDataset();
            obj.createSettings();
        end
        function compute(obj)
            obj.splitTrainingAndTest();
            obj.encodeTensors();
            obj.train();
        end
        function plotFit(obj,m1,m2,s,r)
            d           = load(obj.featureFile);
            isComponent = ~ismember(d.columns,[obj.featureNames,"a","r","rho"]);
            shape       = [d.data(:,ismember(d.columns,["m1","m2","s"])),d.data(:,d.columns == "r")];
            [~,i]       = min(vecnorm(shape-[m1,m2,s,r],2,2));
            rows        = all(abs(shape-shape(i,:)) < 1e-8,2);
            aData       = d.data(rows,d.columns == "a");
            names       = d.columns(isComponent);
            angle       = linspace(-44,134,400)';
            a           = sind(angle)./(sind(angle)+cosd(angle));
            X           = [repmat(shape(i,1:2),numel(a),1),a,repmat(shape(i,3:4),numel(a),1)];
            Y           = obj.loadSurrogate().computeOutputValues(X);
            figure('Color','w')
            for k = 1:numel(names)
                subplot(2,3,k)
                plot(angle,Y(:,k),'-','LineWidth',1.5,'DisplayName','network')
                hold on
                plot(atan2d(aData,1-aData),d.data(rows,find(isComponent,1)+k-1),'k.','MarkerSize',10,'DisplayName','homogenization')
                grid on
                xlabel('wall direction given by a (degrees)')
                title(names(k))
            end
            subplot(2,3,1)
            legend('Location','best')
            sgtitle(sprintf('Network against homogenization, m_1 = %.3g, m_2 = %.3g, s = %.3g, r = %.3g',shape(i,:)))
        end
        function file = getOutputFile(obj)
            file = obj.outputFile;
        end
    end
    methods (Access = private)
        function init(obj)
            obj.datasetFile    = fullfile('latticeParallelogramDataset','dataset.mat');
            obj.testFraction   = 0.15;
            obj.seed           = 1;
            obj.hiddenLayers   = [128 128 128];
            obj.learningRate   = 1e-3;
            obj.maxEpochs      = 8000;
            obj.earlyStop      = 400;
            obj.prefix         = 'latticeParallelogram';
            obj.outputFile     = sprintf('latticeParallelogramNN_Adam_%dx%d_cholesky.mat',obj.hiddenLayers(1),numel(obj.hiddenLayers));
            obj.featureFile    = [obj.prefix,'Features.mat'];
            obj.featureNames   = ["m1","m2","cosA","sinA","s","logR"];
            obj.parameterNames = {'m1','m2','a','s','r'};
        end
        function createFeatureDataset(obj)
            d           = load(obj.datasetFile);
            names       = string(obj.parameterNames);
            params      = zeros(size(d.data,1),numel(names));
            for i = 1:numel(names)
                params(:,i) = d.data(:,d.columns == names(i));
            end
            isComponent  = ~ismember(d.columns,[names,"rho"]);
            f.data       = [ParallelogramCellSurrogate.mapInput(params),params(:,3),params(:,5),d.data(:,d.columns == "rho"),d.data(:,isComponent)];
            f.columns    = [obj.featureNames,"a","r","rho",d.columns(isComponent)];
            f.rhoNominal = d.rhoNominal;
            f.xi         = d.xi;
            f.tensors    = d.tensors;
            f.source     = d.source;
            f.settings   = d.settings;
            save(obj.featureFile,'-struct','f');
        end
        function createSettings(obj)
            d                              = load(obj.featureFile,'columns');
            nFeatures                      = numel(obj.featureNames);
            s                              = latticeSidesNNSettings();
            s.seed                         = obj.seed;
            s.networkParams.hiddenLayers   = obj.hiddenLayers;
            s.optimizerParams.learningRate = obj.learningRate;
            s.optimizerParams.maxEpochs    = obj.maxEpochs;
            s.optimizerParams.earlyStop    = obj.earlyStop;
            s.reference.trainFileName      = [obj.prefix,'Train.mat'];
            s.reference.testFileName       = [obj.prefix,'Test.mat'];
            s.data.trainFileName           = [obj.prefix,'TrainCholesky.mat'];
            s.data.testFileName            = [obj.prefix,'TestCholesky.mat'];
            s.data.xFeatures               = find(ismember(d.columns,obj.featureNames));
            s.data.yFeatures               = find(~ismember(d.columns,[obj.featureNames,"a","r","rho"]));
            s.output.densityFactor         = struct( ...
                'fun',      @(X) 1 - X(:,1).*X(:,2), ...
                'gradient', @(X) [-X(:,2),-X(:,1),zeros(size(X,1),nFeatures-2)]);
            s.transformation               = TurnBoostStretchTransformation();
            s.symmetryColumns              = 1:nFeatures;
            s.parameterNames               = cellstr(obj.featureNames);
            s.evaluation.bins = [ ...
                struct('name','m1','map',@(X) X(:,1),'edges',linspace(0,1,6)), ...
                struct('name','m2','map',@(X) X(:,2),'edges',linspace(0,1,6)), ...
                struct('name','s','map',@(X) X(:,5),'edges',[-Inf,-0.45,-0.15,0.15,0.45,Inf]), ...
                struct('name','log r','map',@(X) X(:,6),'edges',[-Inf,-0.3,-0.1,0.1,0.3,Inf])];
            s.outputFile                   = [obj.prefix,'FeaturesNN.mat'];
            obj.settings                   = s;
        end
        function splitTrainingAndTest(obj)
            splitParallelogramTrainingSet(obj.featureFile,obj.testFraction,obj.seed,obj.prefix);
        end
        function encodeTensors(obj)
            s = obj.settings;
            makeSidesCholeskyDataset(s.reference.trainFileName,s.data.trainFileName,s);
            makeSidesCholeskyDataset(s.reference.testFileName,s.data.testFileName,s);
        end
        function train(obj)
            [problemNN,obj.report,featureSurrogate] = trainSidesNN(obj.settings);
            surrogate      = ParallelogramCellSurrogate(featureSurrogate);
            parameterNames = obj.parameterNames;
            settings       = obj.settings;
            report         = obj.report;
            save(obj.outputFile,'problemNN','surrogate','parameterNames','settings','report');
            obj.surrogate  = surrogate;
            fprintf('Saved %s with inputs (%s)\n',obj.outputFile,strjoin(parameterNames,', '));
        end
        function surrogate = loadSurrogate(obj)
            if isempty(obj.surrogate)
                v             = load(obj.outputFile,'surrogate');
                obj.surrogate = v.surrogate;
            end
            surrogate = obj.surrogate;
        end
    end
end