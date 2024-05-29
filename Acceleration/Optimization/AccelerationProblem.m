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
            s.TOL     = obj.settings.TOL;
            s.maxIter = obj.settings.maxIter;
            s.cost    = obj.cost;
            s.designVariable = obj.designVariable;
            s.momentumParameter.type = 'CONSTANT';
            s.gDescentType = obj.settings.gDescentType;
            t = obj.settings.tau;
            b = obj.settings.beta;
            obj.problem = zeros(numel(t),numel(b));
            Jmin        = obj.problem;
            for i = 1:numel(t)
                for j = 1:numel(b)
                    fprintf('Solving for tau = %.3f and beta = %.2f \n',t(i),b(j))
                    s.tau = t(i);
                    s.momentumParameter.value = b(j);
                    probl = AcceleratedGradientDescent(s);
                    probl.solve();
                    obj.problem(i,j) = numel(probl.J);
                    Jmin(i,j)        = probl.J(end);
                end
            end
            minVal          = min(min(Jmin));
            hasNotConverged = Jmin - minVal > 0.2;
            obj.problem(hasNotConverged) = obj.settings.maxIter;
        end

    end
end