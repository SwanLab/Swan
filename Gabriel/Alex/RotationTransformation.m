classdef RotationTransformation < handle
    methods (Access = public)
        function T = compute(~,params)
            theta = params(1);
            T     = [cosd(theta),-sind(theta);sind(theta),cosd(theta)];
        end
        function Tinv = computeInverse(~,params)
            theta = params(1);
            Tinv  = [cosd(theta),sind(theta);-sind(theta),cosd(theta)];
        end
        function sym = getSymmetry(~,name)
            error('RotationTransformation:UnknownSymmetry','Unknown symmetry: %s',name);
        end
    end
end