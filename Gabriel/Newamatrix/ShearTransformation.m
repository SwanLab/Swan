classdef ShearTransformation < handle

    methods (Access = public)

        function T = compute(~,params)
            a = params(1);
            T = [1,-a;a,1+a^2]/sqrt(1+2*a^2);
        end

        function Tinv = computeInverse(obj,params)
            Tinv = inv(obj.compute(params));
        end

        function sym = getSymmetry(~,name)
            sym.name = name;
            switch name
                case 'mirrorA'
                    sym.mapParams         = @(p) [-p(:,1),p(:,2)];
                    sym.mapParamsJacobian = @(p) repmat([-1 0;0 1],1,1,size(p,1));
                    sym.mapTensor         = @(C) TensorSymmetry.reflect(C,2);
                    sym.isInvariant       = @(p) p(:,1) == 0;
                otherwise
                    error('ShearTransformation:UnknownSymmetry','Unknown symmetry: %s',name);
            end
        end

    end

end