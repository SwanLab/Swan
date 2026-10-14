classdef CostNN < handle 
    
    properties (Access = public)
        value
        gradient
    end

    properties (SetAccess = private, GetAccess = public)
        isBatchDepleted = false
    end

    properties (Access = private)
        shapeFunctions
        weights
        effectiveWeights
        moveBatch

        lossBalancer
        balanceMask
    end

    properties (Access = private)
        shapeValues
    end

    methods (Access = public)

        function obj = CostNN(cParams)

            obj.init(cParams);

        end

        function computeFunctionAndGradient(obj,x)

            compFunc = ...
                @(shI,x) ...
                shI.computeFunctionAndGradient(x);

            [jV,djV] = ...
                obj.computeValueAndGradient(...
                x,compFunc,false);

            obj.value = jV;

            obj.gradient = djV;

        end

        function computeStochasticFunctionAndGradient(obj,x)

            compFunc = ...
                @(shI,x,moveBatch) ...
                shI.computeStochasticCostAndGradient(...
                x,moveBatch);

            [jV,djV] = ...
                obj.computeValueAndGradient(...
                x,compFunc,true);

            obj.value = jV;

            obj.gradient = djV;

        end

        function nF = obtainNumberFields(obj)

            nF = ...
                length(obj.shapeFunctions);

        end

        function titles = getTitleFields(obj)

            nF = ...
                length(obj.shapeFunctions);

            titles = ...
                cell(nF,1);

            for iF = 1:nF

                wI = ...
                    obj.effectiveWeights(iF);

                titleF = ...
                    obj.shapeFunctions{iF}.getTitleToPlot();

                titles{iF} = ...
                    [titleF,...
                     ' (w=',...
                     num2str(wI),...
                     ')'];

            end

        end

        function j = getFields(obj,i)

            j = ...
                obj.shapeValues{i};

        end

        function weights = getEffectiveWeights(obj)

            weights = ...
                obj.effectiveWeights;

        end

        function weights = getBaseWeights(obj)

            weights = ...
                obj.weights;

        end

        function setBatchMover(obj,moveBatch)

            obj.moveBatch = ...
                moveBatch;

        end

        function [alarm,minTestError] = ...
                validateES(obj,alarm,minTestError)

            testError = ...
                obj.shapeFunctions{1}.getTestError();

            if testError < minTestError

                minTestError = ...
                    testError;

                alarm = 0;

            elseif testError == minTestError

                alarm = ...
                    alarm + 0.5;

            else

                alarm = ...
                    alarm + 1;

            end

        end

    end
    
    methods (Access = private)

        function init(obj,cParams)

            obj.shapeFunctions = ...
                cParams.shapeFunctions;

            obj.weights = ...
                cParams.weights(:);

            nF = ...
                length(obj.shapeFunctions);

            if length(obj.weights) ~= nF

                error(...
                    ['Number of CostNN weights must ',...
                     'match number of shape functions.']);

            end

            obj.effectiveWeights = ...
                obj.weights;

            obj.lossBalancer = [];

            obj.balanceMask = ...
                false(nF,1);

            if isfield(cParams,'lossBalancer')

                if ~isempty(cParams.lossBalancer)

                    obj.lossBalancer = ...
                        cParams.lossBalancer;

                    if ~isfield(cParams,'balanceMask')

                        error(...
                            ['balanceMask is required ',...
                             'when a lossBalancer is used.']);

                    end

                    obj.balanceMask = ...
                        logical(...
                        cParams.balanceMask(:));

                    if length(obj.balanceMask) ~= nF

                        error(...
                            ['balanceMask must have ',...
                             'one entry per shape function.']);

                    end

                    if ~any(obj.balanceMask)

                        error(...
                            ['At least one shape function ',...
                             'must be selected for balancing.']);

                    end

                end

            end

        end

        function [jV,djV] = ...
                computeValueAndGradient(...
                obj,x,compFunc,isStochastic)

            nF = ...
                length(obj.shapeFunctions);

            bDa = ...
                false(nF,1);

            Jc = ...
                cell(nF,1);

            dJc = ...
                cell(nF,1);

            for iF = 1:nF

                shI = ...
                    obj.shapeFunctions{iF};

                if isStochastic

                    [j,dJ,bD] = ...
                        compFunc(...
                        shI,x,obj.moveBatch);

                    bDa(iF) = ...
                        bD;

                else

                    [j,dJ] = ...
                        compFunc(...
                        shI,x);

                end

                Jc{iF} = ...
                    j;

                dJc{iF} = ...
                    obj.mergeGradient(dJ);

            end

            obj.updateEffectiveWeights(...
                Jc,isStochastic);

            jV = 0;

            djV = ...
                zeros(size(dJc{1}));

            for iF = 1:nF

                wI = ...
                    obj.effectiveWeights(iF);

                jV = ...
                    jV + wI*Jc{iF};

                djV = ...
                    djV + wI*dJc{iF};

            end

            obj.isBatchDepleted = ...
                any(bDa);

            obj.shapeValues = ...
                Jc;

        end

        function updateEffectiveWeights(...
                obj,Jc,isStochastic)

            obj.effectiveWeights = ...
                obj.weights;

            if isempty(obj.lossBalancer)

                return

            end

            if isStochastic && obj.moveBatch

                lossValues = ...
                    cell2mat(...
                    Jc(obj.balanceMask));

                balancingWeights = ...
                    obj.lossBalancer.update(...
                    lossValues);

            else

                balancingWeights = ...
                    obj.lossBalancer.getWeights();

                if isempty(balancingWeights)

                    balancingWeights = ...
                        ones(...
                        sum(obj.balanceMask),1);

                end

            end

            obj.effectiveWeights(...
                obj.balanceMask) = ...
                obj.weights(...
                obj.balanceMask).*...
                balancingWeights(:);

        end

    end

    methods (Static,Access = private)

        function dJm = mergeGradient(dJ)

            if iscell(dJ)

                nDV = ...
                    length(dJ);

                nDim1 = ...
                    length(dJ{1}.fValues);

                dJm = ...
                    zeros(nDV*nDim1,1);

                for i = 1:nDV

                    ind1 = ...
                        1+nDim1*(i-1);

                    ind2 = ...
                        nDim1+nDim1*(i-1);

                    indices = ...
                        ind1:ind2;

                    dJm(indices,1) = ...
                        dJ{i}.fValues;

                end

            elseif isnumeric(dJ)

                dJm = ...
                    dJ;

            else

                warning(...
                    ['Unsupported input type. ',...
                    'dJ should be a cell array ',...
                    'of structs or a numeric array.']);

            end

        end

    end

end