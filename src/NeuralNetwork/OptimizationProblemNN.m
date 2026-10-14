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
        sensitivity
        regularization
    end
  
    methods (Access = public)

        function obj = OptimizationProblemNN(cParams)

            obj.init(cParams);

            obj.createNetwork();

            obj.createDesignVariable();

            obj.createLossFunctional();

            obj.createSensitivityFunctional();

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
                obj.loss.computeFunctionAndGradient(theta);

            diagnostics.JC = ...
                JC;

            diagnostics.normGC = ...
                norm(gC);

            if obj.useSensitivity

                [JS,gS] = ...
                    obj.sensitivity.computeFunctionAndGradient(...
                    theta);

                diagnostics.JS = ...
                    JS;

                diagnostics.normGS = ...
                    norm(gS);

                denominator = ...
                    max(...
                    norm(gC)*norm(gS),...
                    eps);

                diagnostics.gradientCosine = ...
                    dot(gC,gS)/denominator;

            else

                diagnostics.JS = 0;

                diagnostics.normGS = 0;

                diagnostics.gradientCosine = 0;

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

            fprintf('\nJ_C              = %.12e\n',...
                d.JC)

            fprintf('||grad J_C||     = %.12e\n',...
                d.normGC)

            if obj.useSensitivity

                fprintf('\nJ_sens           = %.12e\n',...
                    d.JS)

                fprintf('||grad J_sens||  = %.12e\n',...
                    d.normGS)

                fprintf('\nGradient cosine  = %.12e\n',...
                    d.gradientCosine)

            end

            fprintf('\nJ_reg            = %.12e\n',...
                d.JR)

            fprintf('||grad J_reg||   = %.12e\n',...
                d.normGR)

            fprintf('\nJ_total          = %.12e\n',...
                d.JTotal)

            fprintf('||grad J_total|| = %.12e\n',...
                d.normGTotal)

            fprintf('\nBase weights\n')
            disp(d.baseWeights')

            fprintf('Effective weights\n')
            disp(d.effectiveWeights')

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

            if isfield(cParams,'sensitivityParams')

                if ~isempty(cParams.sensitivityParams)

                    obj.sensitivityParams = ...
                        cParams.sensitivityParams;

                    obj.useSensitivity = ...
                        true;

                    if ~isfield(...
                            obj.sensitivityParams,...
                            'weight')

                        obj.sensitivityParams.weight = ...
                            1;

                    end

                end

            end

            obj.useBalancing = ...
                false;

            obj.balancingParams = ...
                [];

            obj.lossBalancer = ...
                [];

            if isfield(cParams,'balancingParams')

                if ~isempty(cParams.balancingParams)

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

        function createSensitivityFunctional(obj)

            if ~obj.useSensitivity

                obj.sensitivity = [];

                return

            end

            s.network = ...
                obj.network;

            s.designVariable = ...
                obj.designVariable;

            s.data = ...
                obj.data;

            s.hA = ...
                obj.sensitivityParams.hA;

            s.hB = ...
                obj.sensitivityParams.hB;

            s.sensitivityFileName = ...
                obj.sensitivityParams.sensitivityFileName;

            obj.sensitivity = ...
                SensitivityLossFunctional(s);

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

                s.shapeFunctions = {...
                    obj.loss,...
                    obj.sensitivity,...
                    obj.regularization};

                s.weights = [...
                    1,...
                    obj.sensitivityParams.weight,...
                    obj.costParams.lambda];

            else

                s.shapeFunctions = {...
                    obj.loss,...
                    obj.regularization};

                s.weights = [...
                    1,...
                    obj.costParams.lambda];

            end

            if obj.useBalancing

                if ~isfield(...
                        obj.balancingParams,...
                        'balanceMask')

                    error(...
                        ['balanceMask must be defined ',...
                         'when loss balancing is enabled.']);

                end

                s.lossBalancer = ...
                    obj.lossBalancer;

                s.balanceMask = ...
                    obj.balancingParams.balanceMask;

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
         
    end

end