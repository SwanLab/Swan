classdef QuarterTurnTransformation < handle

    methods (Access = public)

        function T = compute(~,params)
            a = params(1);
            T = [1-a,-a;a,1-a]/sqrt((1-a)^2+a^2);
        end

        function Tinv = computeInverse(~,params)
            a    = params(1);
            Tinv = [1-a,a;-a,1-a]/sqrt((1-a)^2+a^2);
        end

        function sym = getSymmetry(~,name)
            sym.name = name;
            switch name
                case 'mirrorA'
                    sym.mapParams         = @(p) [1-p(:,1),1-p(:,2)];
                    sym.mapParamsJacobian = @(p) repmat([-1 0;0 -1],1,1,size(p,1));
                    sym.mapTensor         = @(C) TensorSymmetry.reflect(C,2);
                    sym.isInvariant       = @(p) p(:,1) == 0.5 & p(:,2) == 0.5;
                otherwise
                    error('QuarterTurnTransformation:UnknownSymmetry','Unknown symmetry: %s',name);
            end
        end

    end

end