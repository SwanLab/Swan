classdef TurnBoostStretchTransformation < handle
    methods (Access = public)
        function T = compute(obj,params)
            [a,s,r] = obj.readParameters(params);
            turn    = [1-a,-a;a,1-a]/sqrt((1-a)^2+a^2);
            boost   = [1,s;s,1]/sqrt(1-s^2);
            stretch = [r,0;0,1/r];
            T       = turn*boost*stretch;
        end
        function Tinv = computeInverse(obj,params)
            Tinv = obj.compute(params)\eye(2);
        end
        function twin = computeTwin(~,params)
            m1   = params(:,1);
            m2   = params(:,2);
            a    = params(:,3);
            s    = params(:,4);
            r    = params(:,5);
            twin = [m2,m1,(1-a)./(1-2*a),-s,1./r];
        end
        function sym = getSymmetry(~,name)
            error('TurnBoostStretchTransformation:UnknownSymmetry','Unknown symmetry: %s',name);
        end
    end
    methods (Static, Access = private)
        function [a,s,r] = readParameters(params)
            a = params(1);
            s = 0;
            r = 1;
            if numel(params) >= 3
                s = params(2);
                r = params(3);
            end
        end
    end
end