classdef Q_learning < handle

    properties (Access = private)
        Qpmax
        e
    end

    methods (Access = public)

        function obj = Q_learning(nFeatures)
            obj.e = zeros(nFeatures, 1);            
        end

        function w = computeW(obj,done,reward,gamma,w,idx,idx_p,Q,alpha,lambda,next_state,params)
            delta = obj.computeDelta(done,reward, gamma, w,Q,next_state,params);
            obj.e(idx) = 1;
            w = w + alpha * delta * obj.e;
            obj.e = obj.computeDecayTrace(gamma,lambda,w,idx_p,done);
        end        

        function delta = computeDelta(obj,done,reward,gamma,w,Q,next_state,params)
            obj.computeQpMax(w,next_state,params)
            if done
                delta = reward - Q;
            else
                delta = reward + gamma * obj.Qpmax - Q;
            end
        end

        function e = computeDecayTrace(obj,gamma,lambda,w,idx_ap,done)


            % Check if chosen action was greedy
            isGreedy = sum(w(idx_ap)) == obj.Qpmax;

            % Reset or decay trace
            if done || ~isGreedy
                e = zeros(size(obj.e));
            else
                e = gamma * lambda * obj.e;
            end
        end
    end

    methods (Access = private)
    
        function computeQpMax(obj,w,next_state,params)
            Qp_all = zeros(1, params.n_actions);
            for a_i = 1:params.n_actions
                idx_p = getActiveTiles(next_state, a_i, params);
                Qp_all(a_i) = sum(w(idx_p));
            end
            Qp_max = max(Qp_all);
            obj.Qpmax = Qp_max;
        end

    end
end