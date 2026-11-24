classdef CellTransformation < handle

    properties (Access = private)
        type
    end

    methods (Access = public)

        function obj = CellTransformation(cParams)
            obj.type = cParams.type;
        end

        function T = compute(obj,params)
            switch obj.type
                case 'BoostSqueeze'
                    T = obj.computeBoostSqueeze(params(1),params(2));
                otherwise
                    error('CellTransformation:UnknownType','Unknown transformation type: %s',obj.type);
            end
        end

        function sym = getSymmetry(obj,name)
            switch obj.type
                case 'BoostSqueeze'
                    sym = obj.getBoostSqueezeSymmetry(name);
                otherwise
                    error('CellTransformation:UnknownType','Unknown transformation type: %s',obj.type);
            end
        end

    end

    methods (Access = private, Static)

        function T = computeBoostSqueeze(a,b)
            boost   = [1,a;a,1]/sqrt(1-a^2);
            squeeze = diag([b,1/b]);
            T       = boost*squeeze;
        end

        function sym = getBoostSqueezeSymmetry(name)
            sym.name = name;
            switch name
                case 'mirrorA'
                    sym.mapParams         = @(p) [-p(:,1),p(:,2)];
                    sym.mapParamsJacobian = @(p) repmat([-1 0;0 1],1,1,size(p,1));
                    sym.mapTensor         = @(C) TensorSymmetry.reflect(C,2);
                    sym.isInvariant       = @(p) p(:,1) == 0;
                case 'inverseB'
                    sym.mapParams         = @(p) [p(:,1),1./p(:,2)];
                    sym.mapParamsJacobian = @(p) reshape([ones(size(p,1),1),zeros(size(p,1),2),-1./p(:,2).^2]',2,2,[]);
                    sym.mapTensor         = @(C) TensorSymmetry.swapAxes(C);
                    sym.isInvariant       = @(p) p(:,2) == 1;
                otherwise
                    error('CellTransformation:UnknownSymmetry','Unknown symmetry: %s',name);
            end
        end

    end

end