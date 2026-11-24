classdef TutorialLatticeNNSidesQuarterTurn < handle
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
    end
    methods (Access = public)
        function obj = TutorialLatticeNNSidesQuarterTurn()
            obj.init();
            obj.createSettings();
        end
        function compute(obj)
            obj.splitTrainingAndTest();
            obj.encodeTensors();
            obj.train();
        end
        function plotFit(obj,m1,m2)
            d           = load(obj.datasetFile);
            isComponent = ~ismember(d.columns,["m1","m2","a","rho"]);
            sides       = d.data(:,ismember(d.columns,["m1","m2"]));
            [~,i]       = min(vecnorm(sides-[m1,m2],2,2));
            rows        = all(sides == sides(i,:),2);
            aData       = d.data(rows,d.columns == "a");
            names       = d.columns(isComponent);
            a           = linspace(min(aData),max(aData),400)';
            X           = [repmat(sides(i,:),numel(a),1),a];
            Y           = obj.loadSurrogate().computeOutputValues(X);
            figure('Color','w')
            for k = 1:numel(names)
                subplot(2,3,k)
                plot(a,Y(:,k),'-','LineWidth',1.5,'DisplayName','network')
                hold on
                plot(aData,d.data(rows,find(isComponent,1)+k-1),'k.','MarkerSize',10,'DisplayName','homogenization')
                grid on
                xlabel('a')
                title(names(k))
            end
            subplot(2,3,1)
            legend('Location','best')
            sgtitle(sprintf('Network against homogenization, m_1 = %.3g, m_2 = %.3g',sides(i,1),sides(i,2)))
        end
        function file = getOutputFile(obj)
            file = obj.outputFile;
        end
    end
    methods (Access = private)
        function init(obj)
            obj.datasetFile  = fullfile('latticeSidesQuarterTurnDataset','dataset.mat');
            obj.testFraction = 0.15;
            obj.seed         = 1;
            obj.hiddenLayers = [64 64 64];
            obj.learningRate = 1e-3;
            obj.maxEpochs    = 100000;
            obj.earlyStop    = 5000;
            obj.prefix       = 'latticeSidesQuarterTurn';
            obj.outputFile   = 'latticeSidesQuarterTurnNN_Adam_64x3_cholesky.mat';
        end
        function createSettings(obj)
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
            s.transformation               = QuarterTurnTransformation();
            s.evaluation.bins = [ ...
                struct('name','m1','map',@(X) X(:,1),'edges',linspace(0,1,6)), ...
                struct('name','m2','map',@(X) X(:,2),'edges',linspace(0,1,6)), ...
                struct('name','a', 'map',@(X) X(:,3),'edges',[-Inf 0 0.5 1 Inf])];
            s.outputFile                   = obj.outputFile;
            obj.settings                   = s;
        end
        function splitTrainingAndTest(obj)
            splitSidesTrainingSetInterior(obj.datasetFile,obj.testFraction,obj.seed,obj.prefix);
        end
        function encodeTensors(obj)
            s = obj.settings;
            makeSidesCholeskyDataset(s.reference.trainFileName,s.data.trainFileName,s);
            makeSidesCholeskyDataset(s.reference.testFileName,s.data.testFileName,s);
        end
        function train(obj)
            [~,obj.report,obj.surrogate] = trainSidesNN(obj.settings);
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