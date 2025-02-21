classdef Cost < handle 
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
        weights
        Msmooth
    end

    properties (Access = private)
        shapeValues
    end

    methods (Access = public)
        function obj = Cost(cParams)
            obj.init(cParams);
        end

        function computeFunctionAndGradient(obj,x)
            nF  = length(obj.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);
            for iF = 1:nF
                shI     = obj.shapeFunctions{iF};
                % Condicional per veure si és stochastic o no
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iF}  = j;
                dJc{iF} = obj.mergeGradient(dJ);
            end
            obj.shapeValues = Jc;
            jV  = 0;
            djV = zeros(size(dJc{1}));
            for iF = 1:nF
                wI  = obj.weights(iF);
                jV  = jV  + wI*Jc{iF};
                djV = djV + wI*dJc{iF};
            end
            obj.value    = jV;
            %obj.gradient = obj.Msmooth*djV;
            obj.gradient = djV;
        end

        function [j,dj,batchesDepleted] = computeStochasticCostAndGradient(obj,x,moveBatch)
            batchesDepleted = obj.shapeFunctions{1}.handleStochasticBatch(moveBatch);
            obj.computeFunctionAndGradient(x);
            j = obj.value;
            dj = obj.gradient;
        end

        function nF = obtainNumberFields(obj)
            nF = length(obj.shapeFunctions);
        end

        function titles = getTitleFields(obj)
            nF = length(obj.shapeFunctions);
            titles = cell(nF,1);
            for iF = 1:nF
                wI         = obj.weights(iF);
                titleF     = obj.shapeFunctions{iF}.getTitleToPlot();
                titles{iF} = [titleF,' (w=',int2str(wI),')'];
            end
        end

        function j = getFields(obj,i)
            j = obj.shapeValues{i};
        end

        % Moure als funcionals --> isEarlyStop()
        function [alarm,minTestError] = validateES(obj,alarm,minTestError)

            [Xtest, Ytest] = obj.shapeFunctions{1}.getTestData();
            [~,y_pred]     = max(loss.getOutput(Xtest),[],2);
            [~,y_target]   = max(Ytest,[],2);

            testError = mean(y_pred ~= y_target);
            if testError < minTestError
                minTestError = testError;
                alarm = 0;
            elseif testError == minTestError
                alarm = alarm + 0.5;
            else
                alarm = alarm + 1;
            end
        end

    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
            obj.weights        = cParams.weights;
            obj.Msmooth        = cParams.Msmooth;
        end

    end

    methods (Static,Access=private)

        function dJm = mergeGradient(dJ)
            if iscell(dJ)
                nDV   = length(dJ);
                nDim1 = length(dJ{1}.fValues);
                dJm   = zeros(nDV*nDim1,1);
                for i = 1:nDV
                    ind1           = 1+nDim1*(i-1);
                    ind2           = nDim1+nDim1*(i-1);
                    indices        = ind1:ind2;
                    dJm(indices,1) = dJ{i}.fValues;
                end
            elseif isnumeric(dJ)
                dJm = dJ;
            else
                warning('Unsupported input type. dJ should be a cell array of structs or a numeric array.');
            end
        end

    end
end
