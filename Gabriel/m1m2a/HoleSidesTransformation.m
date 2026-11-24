classdef HoleSidesTransformation < handle

    methods (Access = public)

        function T = compute(~,params)
            a = params(end);
            T = [1,a;a,1]/sqrt(1-a^2);
        end

        function sym = getSymmetry(~,name)
            sym.name = name;
            switch name
                case 'mirrorA'
                    sym.mapParams         = @(p) [p(:,1),p(:,2),-p(:,3)];
                    sym.mapParamsJacobian = @(p) repmat(diag([1,1,-1]),1,1,size(p,1));
                    sym.mapTensor         = @(C) TensorSymmetry.reflect(C,2);
                    sym.isInvariant       = @(p) p(:,3) == 0;
                case 'swapSides'
                    sym.mapParams         = @(p) [p(:,2),p(:,1),p(:,3)];
                    sym.mapParamsJacobian = @(p) repmat([0 1 0;1 0 0;0 0 1],1,1,size(p,1));
                    sym.mapTensor         = @(C) TensorSymmetry.swapAxes(C);
                    sym.isInvariant       = @(p) p(:,1) == p(:,2);
                otherwise
                    error('HoleSidesTransformation:UnknownSymmetry','Unknown symmetry: %s',name);
            end
        end

    end

end