classdef ConjugateGradientToleranceCalculator < handle

    properties (Access = public)
        val
        oldVal
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
            disp('rho inc norm: ' + string(incNorm));
            t               = obj.tolStandard(incNorm);
            obj.val         = min(obj.tolMax,max(obj.tolMin,t));
        end

        function computeNullSpaceTolerance(obj,type,varargin)
            switch type
                case 'Standard'
                    disp('- New step... -')
                    charNorm = varargin{1};
                    incNorm  = varargin{1};
                    disp('rho inc norm: ' + string(incNorm))                    
                    if incNorm > 0
                        disp('Characteristic norm: ' + string(charNorm))
                        t       = obj.tolStandard(charNorm);%*(1 + 4*incNorm);
                        % obj.val = min(obj.tolMax(charNorm),max(obj.tolMin(charNorm),t));
                        obj.val = min(obj.tolMax(charNorm),t);
                        obj.restartingVal = obj.tolMin(charNorm);
                        obj.updateOldTol();
                    end
                case 'Decreasing'
                    disp('- Decreasing step... -')
                    obj.val = 0.1*obj.oldVal;
                    % obj.updateOldTol();
                case 'Restarting'
                    disp('- Restarting merit function... ')
                    obj.val = obj.restartingVal;
                    % obj.val = 0;
                otherwise
                    error('Wrong type definition')
            end
        end

        function setMMAToleranceParameters(obj)
            % These parameters might be changed in the future
            % obj.tolStandard = @(norm) interp1([0,0.01,0.05,0.21,1,1e5],...
            %     [1e-12,1e-5,1e-5,1e-1,10,10],norm);
            obj.tolStandard = @(norm) 15*norm;
            obj.tolMax = 8e-1;
            obj.tolMin = 1e-8;
        end

        function setNullSpaceToleranceParameters(obj)
            % obj.tolStandard = @(mG) interp1([0,1e-5,1e-3,1e-1,0.2,1,100],[1e-4,1e-3,1e-2,5e-2,0.2,0.08,0.08],mG);
            % obj.tolMax      = @(mG) interp1([0,1e-5,1e-3,1e-1,0.2,1,100],[1e-3,1e-2,5e-2,0.1,0.5,0.1,0.1],mG);
            % obj.tolMin      = @(mG) interp1([0,1e-5,1e-3,1e-1,0.2,1,100],[5e-5,5e-4,1e-3,1e-3,1e-3,1e-3,1e-3],mG);

            % obj.tolStandard = @(mG) interp1([0,1e-5,1e-3,1e-2,0.1,1,100],[1e-4,1e-3,1e-2,5e-2,0.2,0.1,0.1],mG);
            obj.tolStandard = @(mG) 3*mG;
            obj.tolMax      = @(mG) interp1([0,1e-5,1e-3,1e-2,0.1,1,100],[1e-3,1e-2,5e-2,0.1,0.3,0.3,0.3],mG);
            obj.tolMin      = @(mG) interp1([0,1e-5,1e-3,1e-2,0.1,1,100],[5e-5,5e-4,1e-3,1e-3,1e-3,1e-3,1e-3],mG);  

            obj.val         = obj.tolMax(1); % Initial value
        end
        
        function updateOldTol(obj)
            obj.oldVal = obj.val;
        end

    end

    methods (Static, Access = private)

        function emptyFunc(varargin)

        end

    end
end