classdef Data < handle

    properties (Access = public)
        nFeatures
        nSamples
        nLabels
        Xtrain
        Ytrain
        Xtest
        Ytest
        XtrainPhysical
        YtrainPhysical
        XtestPhysical
        YtestPhysical
        Ntest
    end

    properties (Access = private)
        X
        Y
        polynomialOrder
        data
        fileName
        trainFileName
        testFileName
        useSeparatedData
        testRatio
        xFeatures
        yFeatures
        xMin
        xMax
        xRange
        yScale
    end

    methods (Access = public)

        function obj = Data(cParams)

            obj.init(cParams);

            obj.loadData();

            %obj.buildModel();

            if obj.useSeparatedData

                obj.prepareSeparatedData();

            else

                obj.splitdata();

            end

            obj.nLabels = size(obj.Ytrain,2);

            obj.nFeatures = size(obj.Xtrain,2);

            obj.nSamples = ...
                size(obj.Xtrain,1) + ...
                size(obj.Xtest,1);

        end

        function Xnorm = normalizeInput(obj,X)

            Xnorm = ...
                2*(X-obj.xMin)./obj.xRange - 1;

        end

        function Ynorm = normalizeOutput(obj,Y)

            Ynorm = Y./obj.yScale;

        end

        function Y = denormalizeOutput(obj,Ynorm)

            Y = Ynorm.*obj.yScale;

        end

        function J = denormalizeJacobian(obj,Jnorm)

            xFactor = 2./obj.xRange(:);

            yFactor = obj.yScale(:)';

            scale = xFactor*yFactor;

            J = Jnorm.*scale;

        end

        function plotdata(self,i,j)

            gscatter(...
                self.Xtrain(:,i),...
                self.Xtrain(:,j),...
                self.Ytrain,...
                'rgbcmyk','*');

            xlabel(['X',num2str(i)]);

            ylabel(['X',num2str(j)]);

        end

        function plotCorrRow(obj,k)

            x = obj.data(:,1:end-1);

            nf = size(x,2);

            for i = 1:nf

                nexttile((i-1)*nf+k)

                if i == k

                    histogram(x(:,i))

                else

                    gscatter(...
                        obj.data(:,k),...
                        obj.data(:,i),...
                        obj.Y,...
                        'rgbcmyk','*')

                end

            end

        end

        function plotCorrMatrix(obj)

            x = obj.data(:,1:end-1);

            nf = size(x,2);

            figure

            t = tiledlayout(...
                nf,nf,...
                'TileSpacing','Compact');

            title(t,'Features correlation matrix');

            for i = 1:nf

                obj.plotCorrRow(i);

            end

        end

        function updateHyperparameter(obj,h)

            switch h.type

                case 'testRatio'

                    obj.testRatio = h.value;

                    if ~obj.useSeparatedData

                        obj.splitdata();

                    end

                case 'polyGrade'

                    obj.polynomialOrder = h.value;

                    obj.buildModel();

                    if ~obj.useSeparatedData

                        obj.splitdata();

                    end

            end

        end
        function dYnorm = normalizeDerivative(obj,dY,iFeature)

            dYnorm = ...
                dY.*...
                (obj.xRange(iFeature)/2)./...
                obj.yScale;

        end

    end

    methods (Access = private)

        function init(obj,cParams)

            obj.polynomialOrder = ...
                cParams.polynomialOrder;

            obj.xFeatures = ...
                cParams.xFeatures;

            obj.yFeatures = ...
                cParams.yFeatures;

            obj.fileName = '';

            obj.trainFileName = '';

            obj.testFileName = '';

            obj.testRatio = 0;

            if isfield(cParams,'fileName')

                obj.fileName = ...
                    cParams.fileName;

            end

            if isfield(cParams,'trainFileName')

                obj.trainFileName = ...
                    cParams.trainFileName;

            end

            if isfield(cParams,'testFileName')

                obj.testFileName = ...
                    cParams.testFileName;

            end

            if isfield(cParams,'testRatio')

                obj.testRatio = ...
                    cParams.testRatio;

            end

            obj.useSeparatedData = ...
                ~isempty(obj.trainFileName) && ...
                ~isempty(obj.testFileName);

        end

        function loadData(obj)

            if obj.useSeparatedData

                trainData = ...
                    readmatrix(...
                    fullfile(obj.trainFileName));

                testData = ...
                    readmatrix(...
                    fullfile(obj.testFileName));

                obj.XtrainPhysical = ...
                    trainData(:,obj.xFeatures);

                obj.YtrainPhysical = ...
                    trainData(:,obj.yFeatures);

                obj.XtestPhysical = ...
                    testData(:,obj.xFeatures);

                obj.YtestPhysical = ...
                    testData(:,obj.yFeatures);

                obj.data = ...
                    [trainData;testData];

                obj.X = ...
                    obj.data(:,obj.xFeatures);

                obj.Y = ...
                    obj.data(:,obj.yFeatures);

            else

                f = fullfile(obj.fileName);

                obj.data = readmatrix(f);

                obj.X = ...
                    obj.data(:,obj.xFeatures);

                obj.Y = ...
                    obj.data(:,obj.yFeatures);

            end

        end

        function prepareSeparatedData(obj)

            obj.computeNormalization();

            obj.Xtrain = ...
                obj.normalizeInput(...
                obj.XtrainPhysical);

            obj.Xtest = ...
                obj.normalizeInput(...
                obj.XtestPhysical);

            obj.Ytrain = ...
                obj.normalizeOutput(...
                obj.YtrainPhysical);

            obj.Ytest = ...
                obj.normalizeOutput(...
                obj.YtestPhysical);

            obj.Ntest = ...
                size(obj.XtestPhysical,1);

        end

        function Xful = buildModel(obj)

            x = obj.X;

            d = obj.polynomialOrder;

            x1 = x(:,1);

            x2 = x(:,2);

            cont = 1;

            for g = 1:d

                for a = 0:g

                    Xful(:,cont) = ...
                        x2.^a.*x1.^(g-a);

                    cont = cont + 1;

                end

            end

            obj.X = Xful;

        end

        function exponents = ...
                generateExponents(obj,targetDeg)

            exponents = [];

            currentExponents = ...
                zeros(1,obj.nFeatures);

            initialIndex = 1;

            exponents = ...
                obj.generateExponentsRecursive(...
                targetDeg,...
                initialIndex,...
                currentExponents,...
                exponents);

        end

        function exponents = ...
                generateExponentsRecursive(...
                obj,...
                targetDeg,...
                currentIndex,...
                currentExponents,...
                exponents)

            if currentIndex == obj.nFeatures

                currentExponents(currentIndex) = ...
                    targetDeg;

                exponents = ...
                    [exponents;currentExponents];

            else

                for i = 0:targetDeg

                    currentExponents(currentIndex) = i;

                    exponents = ...
                        obj.generateExponentsRecursive(...
                        targetDeg-i,...
                        currentIndex+1,...
                        currentExponents,...
                        exponents);

                end

            end

        end

        function splitdata(obj)

            nD = size(obj.X,1);

            TP = obj.testRatio;

            r = randperm(nD);

            ntest = round(TP/100*nD);

            ntrain = nD-ntest;

            obj.XtrainPhysical = ...
                obj.X(r(1:ntrain),:);

            obj.XtestPhysical = ...
                obj.X(r(ntrain+1:end),:);

            obj.YtrainPhysical = ...
                obj.Y(r(1:ntrain),:);

            obj.YtestPhysical = ...
                obj.Y(r(ntrain+1:end),:);

            obj.computeNormalization();

            obj.Xtrain = ...
                obj.normalizeInput(...
                obj.XtrainPhysical);

            obj.Xtest = ...
                obj.normalizeInput(...
                obj.XtestPhysical);

            obj.Ytrain = ...
                obj.normalizeOutput(...
                obj.YtrainPhysical);

            obj.Ytest = ...
                obj.normalizeOutput(...
                obj.YtestPhysical);

            obj.Ntest = ntest;

        end

        function computeNormalization(obj)

            obj.xMin = ...
                min(obj.XtrainPhysical,[],1);

            obj.xMax = ...
                max(obj.XtrainPhysical,[],1);

            obj.xRange = ...
                obj.xMax-obj.xMin;

            obj.xRange(...
                obj.xRange == 0) = 1;

            obj.yScale = ...
                max(...
                abs(obj.YtrainPhysical),...
                [],1);

            obj.yScale(...
                obj.yScale == 0) = 1;

        end

    end

end