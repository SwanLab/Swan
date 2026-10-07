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
        moveBatch
        shapeValues
    end


    methods (Access = public)

        function obj = CostNN(cParams)

            obj.init(cParams);

        end


        function computeFunctionAndGradient(obj,x)

            %
            % Full-batch evaluation.
            %

            [jV,djV] = ...
                obj.computeValueAndGradient( ...
                    x,false);

            obj.value    = jV;
            obj.gradient = djV;

        end


        function computeStochasticFunctionAndGradient(obj,x)

            %
            % Mini-batch / stochastic evaluation.
            %

            [jV,djV] = ...
                obj.computeValueAndGradient( ...
                    x,true);

            obj.value    = jV;
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
                    obj.weights(iF);

                titleF = ...
                    obj.shapeFunctions{iF}.getTitleToPlot();

                titles{iF} = ...
                    [titleF, ...
                     ' (w=',num2str(wI),')'];

            end

        end


        function j = getFields(obj,i)

            j = ...
                obj.shapeValues{i};

        end


        function setBatchMover(obj,moveBatch)

            obj.moveBatch = ...
                moveBatch;

        end


        function [alarm,minTestError] = ...
                validateES(obj,alarm,minTestError)

            %
            % Validation is performed using only the
            % data-fitting term.
            %
            % For the surrogate:
            %
            % shapeFunctions{1} = LossFunctional
            %
            % Regularization is intentionally excluded from
            % the validation error.
            %

            testError = ...
                obj.shapeFunctions{1}.getTestError();


            %
            % No validation set, or invalid validation value:
            % do not update early stopping.
            %
            if ~isfinite(testError)
                return
            end


            if testError < minTestError

                %
                % Validation improved.
                %
                minTestError = ...
                    testError;

                alarm = 0;


            elseif testError == minTestError

                %
                % Exactly unchanged.
                %
                alarm = ...
                    alarm + 0.5;


            else

                %
                % Validation became worse.
                %
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
                cParams.weights;

            %
            % Safe default.
            %
            % The SGD may change this through setBatchMover().
            %
            obj.moveBatch = false;


            %
            % Basic consistency check.
            %
            if numel(obj.weights) ~= ...
                    numel(obj.shapeFunctions)

                error( ...
                    'CostNN:Weights', ...
                    ['The number of weights must coincide ', ...
                     'with the number of shape functions.']);

            end

        end


        function [jV,djV] = ...
                computeValueAndGradient(obj,x,isStochastic)

            nF = ...
                length(obj.shapeFunctions);


            %
            % Batch-depletion flag associated with each
            % contribution.
            %
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

                    %
                    % Stochastic interface shared by
                    % LossFunctional and Sh_Func_L2norm.
                    %
                    [j,dJ,bD] = ...
                        shI.computeStochasticCostAndGradient( ...
                            x,obj.moveBatch);

                    bDa(iF) = ...
                        logical(bD);

                else

                    %
                    % Full-batch interface shared by
                    % LossFunctional and Sh_Func_L2norm.
                    %
                    [j,dJ] = ...
                        shI.computeFunctionAndGradient(x);

                    bDa(iF) = false;

                end


                Jc{iF} = ...
                    j;


                %
                % Every gradient contribution must have
                % exactly the same number of entries and
                % orientation as the optimization variable x.
                %
                dJc{iF} = ...
                    obj.mergeGradient(dJ,x);

            end


            %
            % Weighted objective:
            %
            % J(theta) = sum_i w_i J_i(theta)
            %
            jV = 0;

            %
            % The gradient belongs to the same vector space
            % as theta.
            %
            djV = ...
                zeros(size(x));


            for iF = 1:nF

                wI = ...
                    obj.weights(iF);

                jV = ...
                    jV ...
                    + ...
                    wI*Jc{iF};

                djV = ...
                    djV ...
                    + ...
                    wI*dJc{iF};

            end


            %
            % For the current NN problem only LossFunctional
            % owns stochastic batches, while the L2
            % regularization always returns false.
            %
            obj.isBatchDepleted = ...
                any(bDa);


            obj.shapeValues = ...
                Jc;

        end

    end


    methods (Static,Access = private)

        function dJm = mergeGradient(dJ,x)

            if iscell(dJ)

                %
                % Generic concatenation of gradients from
                % several design variables.
                %
                % We do NOT assume that all variables have
                % the same number of degrees of freedom.
                %
                nDV = ...
                    numel(dJ);

                gradientParts = ...
                    cell(nDV,1);


                for i = 1:nDV

                    hasValues = ...
                        (isstruct(dJ{i}) && isfield(dJ{i},'fValues')) ...
                        || ...
                        (isobject(dJ{i}) && isprop(dJ{i},'fValues'));

                    if ~hasValues

                        error( ...
                            'CostNN:GradientFormat', ...
                            ['Cell gradient entry %d does not ', ...
                            'contain fValues.'], ...
                            i);

                    end

                    gradientParts{i} = ...
                        dJ{i}.fValues(:);

                end


                dJm = ...
                    vertcat(gradientParts{:});


            elseif isnumeric(dJ)

                %
                % NN gradients and L2 regularization gradients
                % arrive directly as numeric arrays.
                %
                dJm = dJ;


            else

                error( ...
                    'CostNN:UnsupportedGradient', ...
                    ['Unsupported gradient type. ', ...
                     'dJ must be a cell array or ', ...
                     'a numeric array.']);

            end


            %
            % The gradient and the design variable must
            % contain exactly the same number of degrees
            % of freedom.
            %
            if numel(dJm) ~= numel(x)

                error( ...
                    'CostNN:GradientSize', ...
                    ['Gradient has %d entries, while ', ...
                     'the design variable has %d.'], ...
                    numel(dJm), ...
                    numel(x));

            end


            %
            % Preserve exactly the shape used by theta.
            %
            % Example:
            %
            % x  : 1 x n  -> gradient: 1 x n
            % x  : n x 1  -> gradient: n x 1
            %
            dJm = ...
                reshape(dJm,size(x));

        end

    end

end