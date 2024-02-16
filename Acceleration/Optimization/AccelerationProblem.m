classdef AccelerationProblem < handle

    properties (Access = public)
        problem
    end

    properties (Access = private)
        cost
        designVariable
        settings
        tau
    end

    methods (Access = public)

        function obj = AccelerationProblem(cParams)
            obj.init(cParams);
        end

        function createProblemAndCompute(obj)
            s = obj.settings;
            switch s.problemType
                case 'GENERAL'
                    obj.createGeneralProblemAndCompute();
                case 'TAU_BETA'
                    obj.createTauBetaProblemAndCompute();
                otherwise
                    error('Case not implemented yet.')
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.designVariable = cParams.designVariable;
            obj.settings       = cParams.settings;
        end

        function createGeneralProblemAndCompute(obj)
            obj.problem = {};
            n = 1;
            mPolyak   = {'CONSTANT','CONVEX','CONSTANT'};
            mNesterov = {'CONVEX','CONSTANT'};
            vPolyak   = [0,0,1];
            vNesterov = [0,1];
            s.tau     = obj.settings.tau;
            s.TOL     = obj.settings.TOL;
            s.maxIter = obj.settings.maxIter;
            s.cost    = obj.cost;
            s.designVariable = obj.designVariable;

            for i = 1:numel(mPolyak)                
                s.momentumParameter.type  = mPolyak{i};
                s.momentumParameter.value = vPolyak(i);
                s.gDescentType            = 'Polyak';
                probl = AcceleratedGradientDescent(s);
                probl.solve();
                obj.problem{n} = probl;
                n = n + 1;
            end
            for i = 1:numel(mNesterov)
                s.momentumParameter.type  = mNesterov{i};
                s.momentumParameter.value = vNesterov(i);
                s.gDescentType            = 'Nesterov';
                probl = AcceleratedGradientDescent(s);
                probl.solve();
                obj.problem{n} = probl;
                n = n + 1;
            end
            
            s.constraintCase = 'INEQUALITY';
            % mMMA             = {'CONSTANT','CONVEX','CONSTANT'};
            % vMMA             = [0,0,1];
            mMMA = {'CONSTANT'};
            vMMA = 0;
            for i = 1:numel(mMMA)
                s.momentumParameter.type  = mMMA{i};
                s.momentumParameter.value = vMMA(i);
                probl = AcceleratedMMA(s);
                probl.solve();
                obj.problem{n} = probl;
                n = n + 1;
            end

        end

        function createTauBetaProblemAndCompute(obj)
            
        end

    end
end