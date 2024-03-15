classdef ConjugateGradientToleranceCalculator < handle

    properties (Access = public)
        val
        compute
    end

    properties (Access = private)
        tolMax,tolMin,tolStandard
        restartingVal
    end

    methods (Access = public)

        function obj = ConjugateGradientToleranceCalculator(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.val = 1e-1;
            switch cParams.solver
                case 'CONJUGATE GRADIENT'
                    obj.defineToleranceStrategy(cParams);
                otherwise
                    obj.compute = @obj.emptyFunc;
            end
        end

        function defineToleranceStrategy(obj,cParams)
            switch cParams.optimizer
                case 'MMA'
                    obj.compute = @obj.computeMMATolerance;
                    obj.setMMAToleranceParameters();
                case 'Null Space'
                    obj.compute = @obj.computeNullSpaceTolerance;
                    obj.setNullSpaceToleranceParameters();
                otherwise
                    
            end            
        end

        function computeMMATolerance(obj,incNorm)
            disp('rho inc norm: ', string(incNorm));
            t               = obj.tolStandard(incNorm)*incNorm;
            obj.val         = min(obj.tolMax,max(obj.tolMin,t));
        end

        function computeNullSpaceTolerance(obj,type,varargin)
            switch type
                case 'Standard'
                    disp('- New step... -')
                    meritGradNorm = varargin{1};
                    incNorm       = varargin{2};
                    disp('rho inc norm: ' + string(incNorm))                    
                    if incNorm > 0
                        disp('Merit norm: ' + string(meritGradNorm))
                        t       = obj.tolStandard(meritGradNorm)*(1 + incNorm);
                        obj.val = min(obj.tolMax(meritGradNorm),max(obj.tolMin(meritGradNorm),t));
                        obj.restartingVal = obj.tolMin(meritGradNorm);
                    end
                case 'Decreasing'
                    disp('- Decreasing step... -')
                case 'Restarting'
                    disp('- Restarting merit function... ')
                    obj.val = obj.restartingVal;
                otherwise
                    error('Wrong type definition')
            end
        end

        function setMMAToleranceParameters(obj)
            % These parameters might be changed in the future
            obj.tolStandard = @(norm) interp1([0,0.01,0.05,0.21,1,1e5],...
                [1e-12,1e-5,1e-5,1e-1,10,10],norm);
            obj.tolMax = 8e-1;
            obj.tolMin = 1e-8;
        end

        function setNullSpaceToleranceParameters(obj)
            obj.tolStandard = @(mG) interp1([0,0.2,1,100],[5e-4,1e-2,1e-1,1e-1],mG);
            obj.tolMax      = @(mG) interp1([0,0.2,0.5,1,100],[1e-3,5e-2,8e-2,5e-1,5e-1],mG);
            obj.tolMin      = @(mG) interp1([0,0.2,1,100],[1e-6,1e-3,5e-2,5e-2],mG);
            obj.val         = obj.tolMax(1); % Initial value
        end

    end

    methods (Static, Access = private)

        function emptyFunc(varargin)

        end

    end
end