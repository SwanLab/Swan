classdef Sarsa < handle

    properties
        e
    end

    methods (Access = public)

        function obj = Sarsa(nFeatures)
            obj.e = zeros(nFeatures, 1);
        end

        function w = computeW(obj,done,reward,gamma,w,idx,idx_p,Q,alpha,lambda,next_state,params)
            delta = obj.computeDelta(done,reward, gamma, w,idx_p,Q);
            obj.e(idx) = 1;
            w = w + alpha * delta * obj.e;
            obj.e = gamma * lambda * obj.e;
        end

    end

    methods (Access = private)
        
        function delta = computeDelta(obj,done,reward,gamma,w,idx_p,Q)
            Qp = sum(w(idx_p));
            if done
                delta = reward - Q;
            else
                delta = reward + gamma * Qp - Q;
            end
        end

    end
end