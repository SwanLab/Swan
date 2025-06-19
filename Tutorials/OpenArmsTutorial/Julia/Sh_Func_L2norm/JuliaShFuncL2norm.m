classdef JuliaShFuncL2norm < handle
    properties (Access = private)
        data  % Struct returned from Julia constructor (holds designVariable info)
    end

    methods (Access = public)
        function obj = JuliaShFuncL2norm(params)
            % Call Julia constructor and save returned data
            obj.data = callJuliaClass('Sh_Func_L2norm', 'ShFuncL2norm', params);
        end

        function [j, dj, isBD] = computeStochasticCostAndGradient(obj, x, moveBatch)
            % Optional moveBatch argument
            if nargin < 3
                moveBatch = [];
            end

            % Package input parameters
            params.designVariable = struct('thetavec', obj.data.thetavec);
            params.x = double(x(:));
            params.moveBatch = moveBatch;

            % Call Julia method
            result = callJuliaClass('Sh_Func_L2norm', 'computeStochasticCostAndGradient', params);
            
            % Return outputs
            j    = result.j;
            dj   = result.dj;
            isBD = result.isBD;
        end

        function [j, dj] = computeFunctionAndGradient(obj, x)
            params.designVariable = struct('thetavec', obj.data.thetavec);
            params.x = double(x(:));

            result = callJuliaClass('Sh_Func_L2norm', 'computeFunctionAndGradient', params);

            j  = result.j;
            dj = result.dj;
        end
    end
    methods (Access = private)
        function j = computeCost(obj)
            params.designVariable = struct('thetavec', obj.data.thetavec);
            result = callJuliaClass('Sh_Func_L2norm', 'computeCost', params);
            j = result.j;
        end

        function dj = computeGradient(obj)
            params.designVariable = struct('thetavec', obj.data.thetavec);
            result = callJuliaClass('Sh_Func_L2norm', 'computeGradient', params);
            dj = result.dj;
        end
    end
end
