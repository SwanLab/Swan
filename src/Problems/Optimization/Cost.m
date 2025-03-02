classdef Cost < handle 
    
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
        moveBatch
    end

    properties (Access = private)
        shapeValues
    end

    methods (Access = public)

        function obj = Cost(cParams)
            obj.init(cParams);
        end

        function computeFunctionAndGradient(obj,x)
            compFunc = @(shI, x) shI.computeCostAndGradient(x);
            [jV, djV] = obj.computeValueAndGradient(x, compFunc);
            obj.value    = jV;
            obj.gradient = djV;
        end

        function computeStochasticFunctionAndGradient(obj,x)
            compFunc = @(shI, x, moveBatch) shI.computeStochasticCostAndGradient(x, moveBatch);
            [jV, djV] = obj.computeValueAndGradient(x, compFunc);
            obj.value    = jV;
            obj.gradient = djV;
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

        function setBatchMover(obj, moveBatch)
            obj.moveBatch = moveBatch;
        end

        function [alarm,minTestError] = validateES(obj,alarm,minTestError)
            testError = obj.shapeFunctions{1}.getTestError();
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
        end

        function [jV, djV] = computeValueAndGradient(obj, x, compFunc)
            nF  = length(obj.shapeFunctions);
            bDa  = length(obj.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);
            for iF = 1:nF
                shI = obj.shapeFunctions{iF};
                if nargout(compFunc) == 2
                    [j,dJ]  = compFunc(shI,x);
                    bDa(iF) = false;
                else
                    [j,dJ,bD]  = compFunc(shI,x,obj.moveBatch);
                    bDa(iF) = bD;
                end
                Jc{iF}  = j;
                dJc{iF} = obj.mergeGradient(dJ);
            end
            
            jV  = 0;
            djV = zeros(size(dJc{1}));
            for iF = 1:nF
                wI  = obj.weights(iF);
                jV  = jV  + wI*Jc{iF};
                djV = djV + wI*dJc{iF};
            end

            obj.isBatchDepleted = any(bDa);
            obj.shapeValues = Jc;
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
                warning(['Unsupported input type. dJ should be a ' ...
                    'cell array of structs or a numeric array.']);
            end
        end

    end
end
