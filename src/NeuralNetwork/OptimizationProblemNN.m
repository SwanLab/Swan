classdef OptimizationProblemNN < handle

    properties (Access = private)

        data

        networkParams

        optimizerParams

        costParams

        sensitivityParams

        useSensitivity

        balancingParams

        useBalancing

        lossBalancer

        network

        costFunc

        designVariable

        optimizer

        plotter

        loss

        sensitivities

        regularization

    end

    methods (Access = public)

        function obj = OptimizationProblemNN(cParams)

            obj.init(cParams);

            obj.createNetwork();

            obj.createDesignVariable();

            obj.createLossFunctional();

            obj.createSensitivityFunctionals();

            obj.createRegularizationFunctional();

            obj.createLossBalancer();

            obj.createCost();

            obj.createPlotter();

            obj.createOptimizer();

        end

        function solve(obj)

            obj.optimizer.compute();

        end

        function plotCostFnc(obj)

            obj.optimizer.plotCostFunc();

        end

        function plotBoundary(obj,type)

            obj.plotter.plotBoundary(type);

        end

        function plotConections(obj)

            obj.plotter.plotNetworkStatus();

        end

        function plotConfusionMatrix(obj)

            obj.plotter.drawConfusionMat();

        end

        function plotSurface(obj)

            obj.plotter.drawSurfaceResults();

        end

        function plotImage(obj,row)

            obj.plotter.image(row);

        end

        function E = computeError(obj,X,Y)

            Ypred = ...
                obj.computeOutputValues(X);

            E = ...
                norm(Ypred-Y,'fro')/...
                max(norm(Y,'fro'),eps);

        end

        function yOut = computeOutputValues(obj,X)

            Xnorm = ...
                obj.data.normalizeInput(X);

            Ynorm = ...
                obj.network.computeYOut(Xnorm);

            yOut = ...
                obj.data.denormalizeOutput(Ynorm);

        end

        function dY = computeGradient(obj,X)

            Xnorm = ...
                obj.data.normalizeInput(X);

            dYnorm = ...
                obj.network.networkJacobian(Xnorm);

            dY = ...
                obj.data.denormalizeJacobian(dYnorm);

        end

        function [yOut,dY] = ...
                computeOutputValuesAndGradient(obj,X)

            Xnorm = ...
                obj.data.normalizeInput(X);

            [Ynorm,dYnorm] = ...
                obj.network.computeYOutAndJacobian(Xnorm);

            yOut = ...
                obj.data.denormalizeOutput(Ynorm);

            dY = ...
                obj.data.denormalizeJacobian(dYnorm);

        end

        function diagnostics = ...
                computeLossDiagnostics(obj)

            theta = ...
                obj.designVariable.thetavec;

            [JC,gC] = ...
                obj.loss.computeFunctionAndGradient(...
                theta);

            diagnostics.JC = ...
                JC;

            diagnostics.normGC = ...
                norm(gC);

            diagnostics.sensitivityNames = ...
                {};

            diagnostics.JS = ...
                [];

            diagnostics.normGS = ...
                [];

            diagnostics.gradientCosine = ...
                [];

            if obj.useSensitivity

                nSensitivity = ...
                    length(obj.sensitivities);

                diagnostics.sensitivityNames = ...
                    cell(nSensitivity,1);

                diagnostics.JS = ...
                    zeros(nSensitivity,1);

                diagnostics.normGS = ...
                    zeros(nSensitivity,1);

                diagnostics.gradientCosine = ...
                    zeros(nSensitivity,1);

                for iS = 1:nSensitivity

                    [JS,gS] = ...
                        obj.sensitivities{iS}.computeFunctionAndGradient(...
                        theta);

                    diagnostics.sensitivityNames{iS} = ...
                        obj.sensitivityParams.variables(iS).name;

                    diagnostics.JS(iS) = ...
                        JS;

                    diagnostics.normGS(iS) = ...
                        norm(gS);

                    denominator = ...
                        max(...
                        norm(gC)*norm(gS),...
                        eps);

                    diagnostics.gradientCosine(iS) = ...
                        dot(gC,gS)/...
                        denominator;

                end

            end

            [JR,gR] = ...
                obj.regularization.computeFunctionAndGradient(...
                theta);

            diagnostics.JR = ...
                JR;

            diagnostics.normGR = ...
                norm(gR);

            obj.costFunc.computeFunctionAndGradient(...
                theta);

            diagnostics.JTotal = ...
                obj.costFunc.value;

            diagnostics.normGTotal = ...
                norm(obj.costFunc.gradient);

            diagnostics.baseWeights = ...
                obj.costFunc.getBaseWeights();

            diagnostics.effectiveWeights = ...
                obj.costFunc.getEffectiveWeights();

        end

        function printLossDiagnostics(obj)

            d = ...
                obj.computeLossDiagnostics();

            fprintf('\n')
            fprintf('============================================\n')
            fprintf('          NN LOSS DIAGNOSTICS\n')
            fprintf('============================================\n')

            fprintf(...
                '\nJ_C              = %.12e\n',...
                d.JC)

            fprintf(...
                '||grad J_C||     = %.12e\n',...
                d.normGC)

            for iS = 1:length(d.JS)

                name = ...
                    d.sensitivityNames{iS};

                fprintf(...
                    '\nJ_%s              = %.12e\n',...
                    name,...
                    d.JS(iS))

                fprintf(...
                    '||grad J_%s||     = %.12e\n',...
                    name,...
                    d.normGS(iS))

                fprintf(...
                    'cos(grad JC,J%s) = %.12e\n',...
                    name,...
                    d.gradientCosine(iS))

            end

            fprintf(...
                '\nJ_reg            = %.12e\n',...
                d.JR)

            fprintf(...
                '||grad J_reg||   = %.12e\n',...
                d.normGR)

            fprintf(...
                '\nJ_total          = %.12e\n',...
                d.JTotal)

            fprintf(...
                '||grad J_total|| = %.12e\n',...
                d.normGTotal)

            fprintf('\nBase weights\n')

            disp(...
                d.baseWeights')

            fprintf('Effective weights\n')

            disp(...
                d.effectiveWeights')

        end

        function weights = getEffectiveWeights(obj)

            weights = ...
                obj.costFunc.getEffectiveWeights();

        end

        function history = getBalancingHistory(obj)

            history.loss = [];

            history.weights = [];

            history.rho = [];

            if obj.useBalancing

                history.loss = ...
                    obj.lossBalancer.lossHistory;

                history.weights = ...
                    obj.lossBalancer.weightHistory;

                history.rho = ...
                    obj.lossBalancer.rhoHistory;

            end

        end

        function resetLossBalancer(obj)

            if obj.useBalancing

                obj.lossBalancer.reset();

            end

        end

        function [f,gnorm,isBD] = ...
                computeStochasticLossDiagnostics(...
                obj,moveBatch)

            theta = ...
                obj.designVariable.thetavec;

            obj.costFunc.setBatchMover(...
                moveBatch);

            obj.costFunc.computeStochasticFunctionAndGradient(...
                theta);

            f = ...
                obj.costFunc.value;

            gnorm = ...
                norm(obj.costFunc.gradient);

            isBD = ...
                obj.costFunc.isBatchDepleted;

        end

        function diagnostics = ...
                computeSensitivityDiagnostics(obj)

            if ~obj.useSensitivity

                error(...
                    'Sensitivity functional is not enabled.')

            end

            theta = ...
                obj.designVariable.thetavec;

            [~,gC] = ...
                obj.loss.computeFunctionAndGradient(...
                theta);

            nSensitivity = ...
                length(obj.sensitivities);

            diagnostics = ...
                cell(nSensitivity,1);

            for iS = 1:nSensitivity

                d = ...
                    obj.sensitivities{iS}.computeSensitivityDiagnostics(...
                    theta);

                d.normGC = ...
                    norm(gC);

                denominator = ...
                    max(...
                    norm(gC)*norm(d.g),...
                    eps);

                d.gradientCosineC = ...
                    dot(gC,d.g)/...
                    denominator;

                diagnostics{iS} = ...
                    d;

            end

        end

        function printSensitivityDiagnostics(obj)

            diagnostics = ...
                obj.computeSensitivityDiagnostics();

            fprintf('\n')
            fprintf('============================================\n')
            fprintf('      SENSITIVITY LOSS DIAGNOSTICS\n')
            fprintf('============================================\n')

            for iS = 1:length(diagnostics)

                d = ...
                    diagnostics{iS};

                fprintf(...
                    '\nVariable: %s\n',...
                    d.name)

                fprintf(...
                    'Feature index       = %d\n',...
                    d.featureIndex)

                fprintf(...
                    'J                   = %.12e\n',...
                    d.j)

                fprintf(...
                    '||grad J||          = %.12e\n',...
                    d.normG)

                fprintf(...
                    'FD vs analytical    = %.12e\n',...
                    d.fdAnalyticError)

                fprintf(...
                    'FD vs target        = %.12e\n',...
                    d.fdTargetError)

                fprintf(...
                    'Analytic vs target  = %.12e\n',...
                    d.analyticTargetError)

                fprintf(...
                    'cos(grad JC,J)      = %.12e\n',...
                    d.gradientCosineC)

            end

        end

    end

    methods (Access = private)

        function init(obj,cParams)

            obj.data = ...
                cParams.data;

            obj.networkParams = ...
                cParams.networkParams;

            obj.optimizerParams = ...
                cParams.optimizerParams;

            obj.costParams = ...
                cParams.costParams;

            obj.useSensitivity = ...
                false;

            obj.sensitivityParams = ...
                [];

            if isfield(...
                    cParams,...
                    'sensitivityParams')

                if ~isempty(...
                        cParams.sensitivityParams)

                    obj.sensitivityParams = ...
                        cParams.sensitivityParams;

                    if ~isfield(...
                            obj.sensitivityParams,...
                            'sensitivityFileName')

                        error(...
                            ['sensitivityFileName must be ',...
                             'defined in sensitivityParams.']);

                    end

                    if ~isfield(...
                            obj.sensitivityParams,...
                            'variables')

                        error(...
                            ['variables must be defined ',...
                             'in sensitivityParams.']);

                    end

                    if isempty(...
                            obj.sensitivityParams.variables)

                        error(...
                            ['At least one sensitivity ',...
                             'variable must be defined.']);

                    end

                    variables = ...
                        obj.sensitivityParams.variables;

                    nSensitivity = ...
                        length(variables);

                    requiredFields = {...
                        'name',...
                        'featureIndex',...
                        'targetColumns',...
                        'h'};

                    nFeatures = ...
                        size(...
                        obj.data.Xtrain,...
                        2);

                    for iS = 1:nSensitivity

                        for iField = ...
                                1:length(requiredFields)

                            fieldName = ...
                                requiredFields{iField};

                            if ~isfield(...
                                    variables(iS),...
                                    fieldName)

                                error(...
                                    ['Missing field ',...
                                     fieldName,...
                                     ' in sensitivity variable.']);

                            end

                        end

                        featureIndex = ...
                            variables(iS).featureIndex;

                        if ~isscalar(featureIndex) || ...
                                featureIndex < 1 || ...
                                featureIndex > nFeatures || ...
                                featureIndex ~= floor(featureIndex)

                            error(...
                                ['Invalid featureIndex in ',...
                                 'sensitivity variable.']);

                        end

                        if ~isscalar(variables(iS).h) || ...
                                variables(iS).h <= 0

                            error(...
                                ['Sensitivity step h must ',...
                                 'be a positive scalar.']);

                        end

                        if ~isfield(...
                                variables(iS),...
                                'weight')

                            variables(iS).weight = ...
                                1;

                        end

                    end

                    obj.sensitivityParams.variables = ...
                        variables;

                    obj.useSensitivity = ...
                        true;

                end

            end

            obj.useBalancing = ...
                false;

            obj.balancingParams = ...
                [];

            obj.lossBalancer = ...
                [];

            if isfield(...
                    cParams,...
                    'balancingParams')

                if ~isempty(...
                        cParams.balancingParams)

                    obj.balancingParams = ...
                        cParams.balancingParams;

                    if isfield(...
                            obj.balancingParams,...
                            'enabled')

                        obj.useBalancing = ...
                            obj.balancingParams.enabled;

                    else

                        obj.useBalancing = ...
                            true;

                    end

                end

            end

        end

        function createNetwork(obj)

            s = ...
                obj.networkParams;

            s.data = ...
                obj.data;

            obj.network = ...
                Network(s);

        end

        function createDesignVariable(obj)

            obj.designVariable = ...
                obj.network.getLearnableVariables();

        end

        function createLossFunctional(obj)

            s.network = ...
                obj.network;

            s.designVariable = ...
                obj.designVariable;

            s.data = ...
                obj.data;

            s.costType = ...
                obj.costParams.costType;

            obj.loss = ...
                LossFunctional(s);

        end

        function createSensitivityFunctionals(obj)

            if ~obj.useSensitivity

                obj.sensitivities = ...
                    {};

                return

            end

            variables = ...
                obj.sensitivityParams.variables;

            nSensitivity = ...
                length(variables);

            obj.sensitivities = ...
                cell(nSensitivity,1);

            for iS = 1:nSensitivity

                s = ...
                    struct();

                s.network = ...
                    obj.network;

                s.designVariable = ...
                    obj.designVariable;

                s.data = ...
                    obj.data;

                s.sensitivityFileName = ...
                    obj.sensitivityParams.sensitivityFileName;

                s.variableName = ...
                    variables(iS).name;

                s.featureIndex = ...
                    variables(iS).featureIndex;

                s.targetColumns = ...
                    variables(iS).targetColumns;

                s.h = ...
                    variables(iS).h;

                obj.sensitivities{iS} = ...
                    SensitivityLossFunctional(s);

            end

        end

        function createRegularizationFunctional(obj)

            s.designVariable = ...
                obj.designVariable;

            obj.regularization = ...
                Sh_Func_L2norm(s);

        end

        function createLossBalancer(obj)

            if ~obj.useBalancing

                obj.lossBalancer = [];

                return

            end

            if isfield(...
                    obj.balancingParams,...
                    'type')

                type = ...
                    obj.balancingParams.type;

            else

                type = ...
                    'ReLoBRaLo';

            end

            switch lower(type)

                case 'relobralo'

                    obj.lossBalancer = ...
                        ReLoBRaLo(...
                        obj.balancingParams);

                otherwise

                    error(...
                        ['Unknown loss balancing ',...
                         'method: ',type]);

            end

        end

        function createCost(obj)

            if obj.useSensitivity

                nSensitivity = ...
                    length(obj.sensitivities);

                nFunctions = ...
                    nSensitivity + 2;

                s.shapeFunctions = ...
                    cell(nFunctions,1);

                s.weights = ...
                    zeros(nFunctions,1);

                s.shapeFunctions{1} = ...
                    obj.loss;

                s.weights(1) = ...
                    1;

                for iS = 1:nSensitivity

                    index = ...
                        iS + 1;

                    s.shapeFunctions{index} = ...
                        obj.sensitivities{iS};

                    s.weights(index) = ...
                        obj.sensitivityParams.variables(iS).weight;

                end

                s.shapeFunctions{end} = ...
                    obj.regularization;

                s.weights(end) = ...
                    obj.costParams.lambda;

            else

                s.shapeFunctions = {...
                    obj.loss;...
                    obj.regularization};

                s.weights = [...
                    1;...
                    obj.costParams.lambda];

            end

            if obj.useBalancing

                s.lossBalancer = ...
                    obj.lossBalancer;

                if isfield(...
                        obj.balancingParams,...
                        'balanceMask') && ...
                        ~isempty(...
                        obj.balancingParams.balanceMask)

                    s.balanceMask = ...
                        obj.balancingParams.balanceMask;

                else

                    s.balanceMask = ...
                        obj.createDefaultBalanceMask(...
                        length(s.shapeFunctions));

                end

            end

            obj.costFunc = ...
                CostNN(s);

        end

        function createOptimizer(obj)

            s = ...
                obj.optimizerParams;

            s.costFunc = ...
                obj.costFunc;

            s.designVariable = ...
                obj.network.getLearnableVariables();

            s.type = ...
                'SGD';

            s.data = ...
                obj.data;

            s.Xtrain = ...
                obj.data.Xtrain;

            s.Ytrain = ...
                obj.data.Ytrain;

            s.Xtest = ...
                obj.data.Xtest;

            s.Ytest = ...
                obj.data.Ytest;

            s.plotter = ...
                obj.plotter;

            obj.optimizer = ...
                Trainer.create(s);

        end

        function createPlotter(obj)

            s.network = ...
                obj.network;

            s.data = ...
                obj.data;

            s.costfunc = ...
                obj.costFunc;

            obj.plotter = ...
                PlotterNN(s);

        end

        function balanceMask = ...
                createDefaultBalanceMask(obj,nFunctions)

            balanceMask = ...
                true(nFunctions,1);

            balanceMask(end) = ...
                false;

        end

    end

end