classdef TensorSymmetry

    methods (Access = public, Static)

        function C = reflect(C,axis)
            n       = size(C,1);
            s       = ones(n,1);
            s(axis) = -1;
            [i,j,k,l] = ndgrid(s,s,s,s);
            C = C.*(i.*j.*k.*l);
        end

        function C = swapAxes(C)
            p = flip(1:size(C,1));
            C = C(p,p,p,p,:);
        end

    end

end