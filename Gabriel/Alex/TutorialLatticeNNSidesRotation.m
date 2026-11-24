classdef TutorialLatticeNNSidesRotation < handle
    properties (Access = public)
        surrogate
        report
    end
    properties (Access = private)
        datasetFile
        harmonic
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
        function obj = TutorialLatticeNNSidesRotation()
            obj.init();
            obj.createFeatureDataset();
            obj.createSettings();
        end
        function compute(obj)
            obj.splitTrainingAndTest();
            obj.encodeTensors();
            obj.train();
        end
        function plotFit(obj,m1,m2)
            d           = load(obj.featureFile);
            isComponent = ~ismember(d.columns,[obj.featureNames,"theta","rho"]);
            sides       = d.data(:,ismember(d.columns,["m1","m2"]));
            [~,i]       = min(vecnorm(sides-[m1,m2],2,2));
            rows        = all(sides == sides(i,:),2);
            thetaData   = d.data(rows,d.columns == "theta");
            names       = d.columns(isComponent);
            theta       = linspace(min(thetaData),min(thetaData)+d.period,400)';
            X           = [repmat(sides(i,:),numel(theta),1),deg2rad(theta)];
            Y           = obj.loadSurrogate().computeOutputValues(X);
            figure('Color','w')
            for k = 1:numel(names)
                subplot(2,3,k)
                plot(theta,Y(:,k),'-','LineWidth',1.5,'DisplayName','network')
                hold on
                plot(thetaData,d.data(rows,find(isComponent,1)+k-1),'k.','MarkerSize',10,'DisplayName','homogenization')
                grid on
                xlabel('\theta (degrees)')
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
            obj.datasetFile    = fullfile('latticeSidesRotationDataset','datasetFullTurn.mat');
            obj.harmonic       = 2;
            obj.testFraction   = 0.15;
            obj.seed           = 1;
            obj.hiddenLayers   = [64 64 64];
            obj.learningRate   = 1e-3;
            obj.maxEpochs      = 100000;
            obj.earlyStop      = 5000;
            obj.prefix         = 'latticeSidesRotation';
            obj.outputFile     = 'latticeSidesRotationNN_Adam_64x3_cholesky.mat';
            obj.featureFile    = [obj.prefix,'Features.mat'];
            obj.featureNames   = ["m1","m2","cosTheta","sinTheta"];
            obj.parameterNames = {'m1','m2','theta'};
        end
        function createFeatureDataset(obj)
            d           = load(obj.datasetFile);
            m1          = d.data(:,d.columns == "m1");
            m2          = d.data(:,d.columns == "m2");
            theta       = d.data(:,d.columns == "theta");
            isComponent = ~ismember(d.columns,["m1","m2","theta","rho"]);
            angle       = obj.harmonic*deg2rad(theta);
            f.data       = [m1,m2,cos(angle),sin(angle),theta,d.data(:,d.columns == "rho"),d.data(:,isComponent)];
            f.columns    = [obj.featureNames,"theta","rho",d.columns(isComponent)];
            f.rhoNominal = 1 - m1.*m2;
            f.xi         = atan2(m1,m2);
            f.tensors    = d.tensors;
            f.source     = repmat("computed",numel(theta),1);
            f.period     = d.period;
            f.settings   = d.settings;
            save(obj.featureFile,'-struct','f');
        end
        function createSettings(obj)
            d                              = load(obj.featureFile,'columns');
            nFeatures                      = numel(obj.featureNames);
            harmonic                       = obj.harmonic;
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
            s.data.yFeatures               = find(~ismember(d.columns,[obj.featureNames,"theta","rho"]));
            s.output.densityFactor         = struct( ...
                'fun',      @(X) 1 - X(:,1).*X(:,2), ...
                'gradient', @(X) [-X(:,2),-X(:,1),zeros(size(X,1),nFeatures-2)]);
            s.transformation               = RotationTransformation();
            s.symmetryColumns              = 1:nFeatures;
            s.parameterNames               = cellstr(obj.featureNames);
            s.evaluation.bins = [ ...
                struct('name','m1','map',@(X) X(:,1),'edges',linspace(0,1,6)), ...
                struct('name','m2','map',@(X) X(:,2),'edges',linspace(0,1,6)), ...
                struct('name','theta (degrees)','map',@(X) atan2d(X(:,4),X(:,3))/harmonic,'edges',linspace(-90,90,5))];
            s.outputFile                   = [obj.prefix,'FeaturesNN.mat'];
            obj.settings                   = s;
        end
        function splitTrainingAndTest(obj)
            splitSidesRotationTrainingSet(obj.featureFile,obj.testFraction,obj.seed,obj.prefix);
        end
        function encodeTensors(obj)
            s = obj.settings;
            makeSidesCholeskyDataset(s.reference.trainFileName,s.data.trainFileName,s);
            makeSidesCholeskyDataset(s.reference.testFileName,s.data.testFileName,s);
        end
        function train(obj)
            [problemNN,obj.report,featureSurrogate] = trainSidesNN(obj.settings);
            surrogate      = PeriodicAngleSurrogate(featureSurrogate,obj.harmonic);
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