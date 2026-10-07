classdef Sh_Func_L2norm < handle

    properties (Access = private)
        designVariable
    end


    methods (Access = public)

        function obj = Sh_Func_L2norm(cParams)

            obj.init(cParams);

        end


        function [j,dj,isBD] = ...
                computeStochasticCostAndGradient(obj,x,~)

            %
            % The L2 regularization is deterministic.
            % It does not depend on the mini-batch.
            %

            [j,dj] = ...
                obj.computeFunctionAndGradient(x);

            isBD = false;

        end
        function t = getTitleToPlot(~)

            t = 'L2 regularization';

        end


        function [j,dj] = ...
                computeFunctionAndGradient(obj,x)

            obj.designVariable.thetavec = x;

            j  = obj.computeCost();
            dj = obj.computeGradient();

        end

    end


    methods (Access = private)

        function init(obj,cParams)

            obj.designVariable = ...
                cParams.designVariable;

        end


        function j = computeCost(obj)

            theta = ...
                obj.designVariable.thetavec;

            %
            % R(theta) = 1/2 ||theta||_2^2
            %
            % Written independently of whether theta is
            % stored as a row or column vector.
            %
            j = ...
                0.5 * sum(theta(:).^2);

        end


        function dj = computeGradient(obj)

            theta = ...
                obj.designVariable.thetavec;

            %
            % dR/dtheta = theta
            %
            % Preserve exactly the orientation of theta.
            %
            dj = theta;

        end

    end

end