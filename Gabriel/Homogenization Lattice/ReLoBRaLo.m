classdef ReLoBRaLo < handle

    properties (SetAccess = private, GetAccess = public)
        weightHistory
        lossHistory
        rhoHistory
    end

    properties (Access = private)
        alpha
        temperature
        rhoExpected
        epsilon

        initialized
        stepCount

        initialLoss
        previousLoss

        weights
    end

    methods (Access = public)

        function obj = ReLoBRaLo(cParams)

            obj.init(cParams);

        end

        function weights = update(obj,lossValues)

            lossValues = lossValues(:);

            if ~obj.initialized

                nTerms = length(lossValues);

                obj.weights = ...
                    ones(nTerms,1);

                obj.previousLoss = ...
                    ones(nTerms,1);

                obj.initialLoss = ...
                    ones(nTerms,1);

                obj.initialized = true;

            end

            if length(lossValues) ~= ...
                    length(obj.weights)

                error(...
                    ['Number of loss terms changed ',...
                     'during ReLoBRaLo training.']);

            end

            lambdaPrevious = ...
                obj.computeRelativeWeights(...
                lossValues,...
                obj.previousLoss);

            lambdaInitial = ...
                obj.computeRelativeWeights(...
                lossValues,...
                obj.initialLoss);

            rho = ...
                double(...
                rand < obj.rhoExpected);

            if obj.stepCount == 0

                alphaCurrent = 1;

            elseif obj.stepCount == 1

                alphaCurrent = 0;

            else

                alphaCurrent = ...
                    obj.alpha;

            end

            obj.weights = ...
                rho*alphaCurrent*obj.weights + ...
                (1-rho)*alphaCurrent*lambdaInitial + ...
                (1-alphaCurrent)*lambdaPrevious;

            obj.previousLoss = ...
                lossValues;

            if obj.stepCount == 1

                obj.initialLoss = ...
                    lossValues;

            end

            obj.storeState(...
                lossValues,rho);

            obj.stepCount = ...
                obj.stepCount + 1;

            weights = ...
                obj.weights;

        end

        function weights = getWeights(obj)

            weights = ...
                obj.weights;

        end

        function reset(obj)

            obj.initialized = false;

            obj.stepCount = 0;

            obj.initialLoss = [];

            obj.previousLoss = [];

            obj.weights = [];

            obj.weightHistory = [];

            obj.lossHistory = [];

            obj.rhoHistory = [];

        end

    end

    methods (Access = private)

        function init(obj,cParams)

            obj.alpha = ...
                cParams.alpha;

            obj.temperature = ...
                cParams.temperature;

            obj.rhoExpected = ...
                cParams.rhoExpected;

            if isfield(cParams,'epsilon')

                obj.epsilon = ...
                    cParams.epsilon;

            else

                obj.epsilon = ...
                    1e-12;

            end

            obj.initialized = false;

            obj.stepCount = 0;

            obj.initialLoss = [];

            obj.previousLoss = [];

            obj.weights = [];

            obj.weightHistory = [];

            obj.lossHistory = [];

            obj.rhoHistory = [];

        end

        function lambda = ...
                computeRelativeWeights(...
                obj,currentLoss,referenceLoss)

            nTerms = ...
                length(currentLoss);

            z = ...
                currentLoss./...
                (obj.temperature*referenceLoss + ...
                obj.epsilon);

            z = ...
                z-max(z);

            exponential = ...
                exp(z);

            lambda = ...
                nTerms*...
                exponential/...
                sum(exponential);

        end

        function storeState(...
                obj,lossValues,rho)

            obj.lossHistory(:,end+1) = ...
                lossValues;

            obj.weightHistory(:,end+1) = ...
                obj.weights;

            obj.rhoHistory(1,end+1) = ...
                rho;

        end

    end

end