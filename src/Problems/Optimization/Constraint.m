classdef Constraint < handle

    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
        Msmooth
    end

    methods (Access = public)

        function obj = Constraint(cParams)
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

            jV     = zeros(nF,1);
            nTotal = numel(dJc{1});
            djV    = zeros(nTotal,nF);

            for iF = 1:nF
                if numel(dJc{iF}) ~= nTotal
                    error('Constraint:InconsistentGradientSize', ...
                        'All constraint gradients must have the same size.');
                end

                jV(iF)    = Jc{iF};
                djV(:,iF) = dJc{iF};
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
                titles{iF} = obj.shapeFunctions{iF}.getTitleToPlot();
            end
        end

    end

    methods (Access = private)

        function obj = init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
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

            dJm    = zeros(nTotal,1);
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