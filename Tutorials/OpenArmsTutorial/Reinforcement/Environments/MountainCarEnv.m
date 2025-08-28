classdef MountainCarEnv
    properties (Access = private)
        posMin
        posMax
        velMin
        velMax
    end
    
    methods (Access = public)
        function obj = MountainCarEnv(params)
            posRange = params.pos_range;
            velRange = params.vel_range;
            obj.posMin = posRange(1);
            obj.posMax = posRange(2);
            obj.velMin = velRange(1);
            obj.velMax = velRange(2);
        end

        function state = reset(obj, s0)
            if nargin == 2
                state = s0;
            else
                pos = obj.posMin + (obj.posMax - obj.posMin) * rand();
                vel = obj.velMin + (obj.velMax - obj.velMin) * rand();
                state = [pos; vel];
            end
        end

        function [nextState, reward, done] = step(obj, state, action)
            pos = state(1);
            vel = state(2);

            % Convert action: 1=left, 2=neutral, 3=right
            force = (action - 2) * 0.001;
            gravity = -0.0025 * cos(3 * pos);

            % Velocity update
            vel = vel + force + gravity;
            vel = min(max(vel, obj.velMin), obj.velMax);

            % Position update
            pos = pos + vel;
            pos = min(max(pos, obj.posMin), obj.posMax);

            if pos == obj.posMin
                vel = 0;
            end

            nextState = [pos; vel];
            done = pos >= 0.5;
            reward = -1;
        end
    end
end
