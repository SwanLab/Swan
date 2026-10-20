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
            nF  = numel(obj.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);

            for iF = 1:nF
                shI     = obj.shapeFunctions{iF};
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iF}  = j;
                dJc{iF} = obj.mergeGradient(dJ);
            end

            obj.shapeValues = Jc;

            jV  = 0;
            djV = zeros(size(dJc{1}));

            for iF = 1:nF
                wI  = obj.weights(iF);
                jV  = jV + wI*Jc{iF};
                djV = djV + wI*dJc{iF};
            end

            obj.value    = jV;
            obj.gradient = obj.Msmooth*djV;
        end

        function nF = obtainNumberFields(obj)
            nF = numel(obj.shapeFunctions);
        end

        function titles = getTitleFields(obj)
            nF     = numel(obj.shapeFunctions);
            titles = cell(nF,1);

            for iF = 1:nF
                wI         = obj.weights(iF);
                titleF     = obj.shapeFunctions{iF}.getTitleToPlot();
                titles{iF} = [titleF,' (w=',num2str(wI),')'];
            end
        end

        function j = getFields(obj,i)
            j = cell2mat(obj.shapeValues(i));
        end

    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
            obj.weights        = cParams.weights;
            obj.Msmooth        = cParams.Msmooth;
        end

    end

    methods (Static, Access = private)

        function dJm = mergeGradient(dJ)
            nDV    = numel(dJ);
            nTotal = 0;

            for i = 1:nDV
                nTotal = nTotal + numel(dJ{i}.fValues);
            end

            dJm   = zeros(nTotal,1);
            offset = 0;

            for i = 1:nDV
                values  = dJ{i}.fValues(:);
                nValues = numel(values);
                indices = offset + (1:nValues);

                dJm(indices) = values;

                offset = offset + nValues;
            end
        end

    end

end